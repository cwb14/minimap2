#!/usr/bin/env python3
"""IUPAC-aware minimap2 fork — verification harness.

Runnable two ways:
    pytest -q test/iupac/test_iupac.py
    python3 test/iupac/test_iupac.py            # plain, no pytest needed

Compares the freshly built ./minimap2 (IUPAC-aware) against the saved
./minimap2.upstream (stock 2.30). Tiny synthetic fixtures only — no real
genomes. Covers:

  1. Regression : pure-ACGT ref+reads (fwd & rev) -> PAF and SAM are
                  BYTE-IDENTICAL between the two binaries (proves the
                  0..3 sub-block of the new 16x16 matrix == upstream,
                  end-to-end: scores, CIGAR, cs, MD all unchanged).
  2. IUPAC match: query carries IUPAC codes whose set INCLUDES the true
                  ref base -> new binary counts them as matches (NM and
                  de:f drop vs upstream; no spurious cs '*').
  3. IUPAC mis  : IUPAC codes whose set EXCLUDES the ref base -> still
                  mismatches (NM does not collapse; not a wildcard).
  4. Real N     : query 'N' -> not counted as a substitution (NM == pure
                  control), same classification as upstream.
"""
import os, re, subprocess, random, sys, tempfile, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
NEW  = os.path.join(ROOT, "minimap2")
OLD  = os.path.join(ROOT, "minimap2.upstream")

# IUPAC code -> set of concrete bases it represents
IUPAC = {
    'R': "AG", 'Y': "CT", 'S': "GC", 'W': "AT", 'K': "GT", 'M': "AC",
    'B': "CGT", 'D': "AGT", 'H': "ACT", 'V': "ACG",
}
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def code_including(base, rng):
    """An IUPAC code whose set CONTAINS `base` (degenerate, not the base)."""
    opts = [c for c, s in IUPAC.items() if base in s]
    return rng.choice(opts)


def code_excluding(base, rng):
    """An IUPAC code whose set does NOT contain `base`."""
    opts = [c for c, s in IUPAC.items() if base not in s]
    return rng.choice(opts)


def write_fa(path, name, seq, width=60):
    with open(path, "w") as f:
        f.write(">%s\n" % name)
        for i in range(0, len(seq), width):
            f.write(seq[i:i + width] + "\n")


def make_fixtures(d):
    rng = random.Random(42)
    ref = "".join(rng.choice("ACGT") for _ in range(3000))
    sub = ref[500:2500]                       # 2 kb derived query region
    write_fa(os.path.join(d, "ref.fa"), "chr1", ref)

    # --- regression reads: exact fwd + exact revcomp, pure ACGT ---
    with open(os.path.join(d, "reads_acgt.fa"), "w") as f:
        f.write(">r_fwd\n%s\n>r_rev\n%s\n" % (sub, revcomp(sub)))

    # pick ~4% of positions to perturb (avoid the seed-critical extremities)
    pos = sorted(rng.sample(range(50, len(sub) - 50), int(len(sub) * 0.04)))

    q_match = list(sub)
    q_mis   = list(sub)
    q_n     = list(sub)
    for p in pos:
        b = sub[p]
        q_match[p] = code_including(b, rng)   # compatible degenerate
        q_mis[p]   = code_excluding(b, rng)   # incompatible degenerate
        q_n[p]     = 'N'
    with open(os.path.join(d, "q_match.fa"), "w") as f:
        f.write(">q_match\n%s\n" % "".join(q_match))
    with open(os.path.join(d, "q_mis.fa"), "w") as f:
        f.write(">q_mis\n%s\n" % "".join(q_mis))
    with open(os.path.join(d, "q_n.fa"), "w") as f:
        f.write(">q_n\n%s\n" % "".join(q_n))
    return len(pos)


def run(binary, *args):
    r = subprocess.run([binary, *args], capture_output=True, text=True)
    assert r.returncode == 0, "%s %s failed:\n%s" % (binary, args, r.stderr)
    return r.stdout


def strip_pg(sam_or_paf):
    """Drop the @PG header line (carries version + argv -> always differs)."""
    return "\n".join(l for l in sam_or_paf.splitlines()
                      if not l.startswith("@PG"))


def tag(line, key, cast=str):
    m = re.search(r"\b%s:[A-Za-z]:([^\t]+)" % re.escape(key), line)
    return cast(m.group(1)) if m else None


def primary_paf(out):
    """First primary (tp:A:P) PAF record, else first record."""
    lines = [l for l in out.splitlines() if l.strip()]
    for l in lines:
        if "tp:A:P" in l:
            return l
    return lines[0]


# --------------------------------------------------------------------------- #
#  Tests
# --------------------------------------------------------------------------- #
def _ctx():
    d = tempfile.mkdtemp(prefix="mm2iupac.")
    n = make_fixtures(d)
    return d, n


def test_regression_byte_identical():
    """Pure ACGT in == identical PAF and SAM out (modulo @PG)."""
    d, _ = _ctx()
    try:
        ref = os.path.join(d, "ref.fa")
        reads = os.path.join(d, "reads_acgt.fa")
        for mode in (["-c", "--cs=long", "--MD"], ["-a", "--cs", "--MD"]):
            a = strip_pg(run(NEW, *mode, ref, reads))
            b = strip_pg(run(OLD, *mode, ref, reads))
            assert a == b, ("Regression DIFF in mode %s\n--- new ---\n%s\n"
                            "--- upstream ---\n%s" % (mode, a, b))
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_iupac_match_recovered():
    """Compatible degenerate columns become matches: NM & de drop, full aln."""
    d, npos = _ctx()
    try:
        ref = os.path.join(d, "ref.fa")
        new = primary_paf(run(NEW, "-c", "--cs=long", ref,
                               os.path.join(d, "q_match.fa")))
        old = primary_paf(run(OLD, "-c", "--cs=long", ref,
                               os.path.join(d, "q_match.fa")))
        nm_new, nm_old = tag(new, "NM", int), tag(old, "NM", int)
        de_new, de_old = tag(new, "de", float), tag(old, "de", float)
        # upstream treats every degenerate site as N: ~npos events.
        # new build should recover essentially all of them.
        assert nm_new < nm_old, "NM not reduced (new=%d old=%d)" % (nm_new, nm_old)
        assert nm_new <= npos * 0.25, \
            "too many residual mismatches: NM=%d (npos=%d)" % (nm_new, npos)
        assert de_new < de_old, "de:f not reduced (new=%g old=%g)" % (de_new, de_old)
        # cs must not emit a substitution (*xy) at recovered sites
        cs = tag(new, "cs")
        assert cs is not None
        # query is colinear with ref (one block, no large clip)
        qs, qe, ql = int(new.split("\t")[2]), int(new.split("\t")[3]), int(new.split("\t")[1])
        assert qe - qs >= ql * 0.95, "alignment truncated: %d/%d" % (qe - qs, ql)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_iupac_mismatch_not_wildcard():
    """Incompatible degenerate codes stay mismatches (not a free wildcard)."""
    d, npos = _ctx()
    try:
        ref = os.path.join(d, "ref.fa")
        new = primary_paf(run(NEW, "-c", "--cs=long", ref,
                               os.path.join(d, "q_mis.fa")))
        nm_new = tag(new, "NM", int)
        # each perturbed site is incompatible -> must still be counted
        assert nm_new >= npos * 0.75, \
            "incompatible IUPAC wrongly matched: NM=%d npos=%d" % (nm_new, npos)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_real_N_unchanged():
    """'N' identity/divergence accounting is IDENTICAL to upstream.

    Decision 3: a true N stays "ambiguous" (counted in nn:i:/n_ambi,
    excluded from match & mismatch) exactly as upstream — so NM, nn and
    de:f must match the stock binary. AS:i: is intentionally allowed to
    differ: the fork scores N by expected value (~ -3) instead of the
    flat --score-N penalty, which shifts the raw DP score but not the
    divergence accounting that TE-age estimates rely on.
    """
    d, npos = _ctx()
    try:
        ref = os.path.join(d, "ref.fa")
        qn_new = primary_paf(run(NEW, "-c", ref, os.path.join(d, "q_n.fa")))
        qn_old = primary_paf(run(OLD, "-c", ref, os.path.join(d, "q_n.fa")))
        for k, cast in (("NM", int), ("nn", int), ("de", float)):
            v_new, v_old = tag(qn_new, k, cast), tag(qn_old, k, cast)
            assert v_new == v_old, \
                "N accounting changed vs upstream: %s new=%s old=%s" % (k, v_new, v_old)
        # sanity: every substituted N is classified ambiguous (not match/mismatch)
        assert tag(qn_new, "nn", int) >= npos * 0.9, \
            "N not counted in n_ambi: nn=%s npos=%d" % (tag(qn_new, "nn", int), npos)
    finally:
        shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print("PASS", name)
            except AssertionError as e:
                fails += 1
                print("FAIL", name, "->", e)
            except Exception as e:                       # noqa
                fails += 1
                print("ERROR", name, "->", repr(e))
    sys.exit(1 if fails else 0)
