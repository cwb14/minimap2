# memo — IUPAC-aware minimap2 fork

- **Date:** 2026-05-18
- **Author/context:** Claude Code session (resume below). Branch `iupac-aware`.
- **Repo:** `/data2/chris/fungi/best/soloLTR/crori/sim/minimap2`
  (own git repo; forked from upstream minimap2 `2.30-r1299`,
  base commit `e5066c7`). New version string: `2.30-r1299-iupac`.

## Purpose

Make minimap2's **alignment stage** interpret IUPAC degeneracy codes
(`R Y S W K M B D H V` + `N`) instead of collapsing every non-ACGT base to
one penalised sentinel. Motivation: TE/LTR **consensus** sequences carry
IUPAC codes at variable columns; treating them as `N` discards real
information and inflates the divergence used for soloLTR / element-age
estimates. Seeding, chaining, sdust, splice and pure-ACGT behaviour are
left **bit-identical to upstream** (regression-safe).

Locked design decisions (with the user):

1. **Scoring** = expected substitution score over the degeneracy:
   `s(x,y) = mean over i∈set(x), j∈set(y) of base(i,j)`. N (={ACGT})
   stays penalised (~`(a+3b)/4`).
2. **Seeding** = IUPAC-aware (added in a second commit, default-on). Real
   degeneracy codes (R Y S W K M B D H V) are projected to a canonical
   concrete base (lowest set bit, A<C<G<T) in `mm_sketch` so minimizers
   span degenerate columns; N and `?` still break the k-mer run (genuine
   unknown, not informative — never seed across assembly N-gaps). Pure
   ACGT is untouched (0..3 map to themselves) so the byte-identical
   regression guarantee holds. This is what actually recovers TE-copy
   alignments the IUPAC-blind seeding drops.
3. **Identity/NM/de** = a real-IUPAC column counts as a match (NM+0) iff
   the base sets overlap, else mismatch (NM+1); `N`/unknown keep upstream
   accounting (excluded from match/mismatch, counted in `nn:i:`/`n_ambi`).

## Environment

- OS: `Linux 5.15.0-151-generic x86_64`
- Compiler: `cc (Ubuntu 11.4.0-1ubuntu1~22.04.3) 11.4.0`
- Python: `3.10.18`, `pytest 8.3.4`
- No new dependencies. Build = stock `Makefile` (`-O2 -Wall -Wc++-compat`).

## What changed (source)

| File | Change |
|------|--------|
| `sketch.c` | Extended single `seq_nt4_table` to 16 symbols (0-3 = ACGT unchanged; 4=N, 5..14 = R Y S W K M B D H V, 15=?). Added `mm_seq_nt16_set`, `mm_comp_table`, `mm_seq_nt16_str`/`_lc`. **IUPAC-aware seeding**: `mm_seq_nt4_proj[16]` projects real degeneracy → concrete base; `mm_sketch` uses it (N/? still break). |
| `mmpriv.h` | Extern decls + `static inline mm_iupac_compat()`. |
| `align.c` | `ksw_gen_*_mat` → `mm_gen_iupac_mat` (16×16 expected-score). All `mat[25]`→`[256]`, alphabet size `5`→`MM_ASIZE(16)` at every ksw call/profile. IUPAC revcomp via `mm_comp_table`. Identity/`n_ambi` stats: N/? excluded as upstream, real IUPAC = compat?match:mismatch. Ungapped SR path uses the matrix. |
| `index.c` | `mm_idx_getseq_rev` reverse-complement → `mm_comp_table` (IUPAC-correct minus strand). Ref store now keeps 0..15 (4-bit packing unchanged). |
| `format.c` | cs/MD/ds match decision via `mm_fmt_is_sub` (consistent with NM/de); base decode via `mm_seq_nt16_str[_lc]` (fixes 5..15 OOB on the old 5-char `"ACGTN"`); splice/qseq revcomp via `mm_comp_table`. |
| `jump.c` | revcomp via `mm_comp_table` (splice path; consistency). |
| `python/cmappy.h` | mappy seq decode → 16-char literal. |
| `minimap.h` | `MM_VERSION` → `2.30-r1299-iupac`. |

`sdust.c` and `bseq.c` `seq_comp_table` need **no change** (sdust's local
table is `#if`-guarded to the standalone binary only; `seq_comp_table` is
already fully IUPAC-aware at the ASCII level upstream).

## Build

```
cd /data2/chris/fungi/best/soloLTR/crori/sim/minimap2
cp -p minimap2 minimap2.upstream   # one-time: keep stock binary for diff
make clean && make -j4
./minimap2 --version               # -> 2.30-r1299-iupac
```

## Verification (all pass)

```
python3 -m pytest -q test/iupac/test_iupac.py        # 4 passed
# or, without pytest:
python3 test/iupac/test_iupac.py
```

Synthetic (seeded, tiny — `test/iupac/test_iupac.py`):

- **Regression** — pure-ACGT ref+reads (fwd & rev), `-c --cs=long --MD`
  and `-a --cs --MD`: new vs `minimap2.upstream` **byte-identical**.
- **IUPAC match** — query with codes whose set includes the ref base:
  `NM 80→0`, `de 0.04→0.0`, full-length, no spurious cs `*`.
- **IUPAC mismatch** — codes excluding the ref base: still mismatches
  (not a wildcard).
- **Real N** — `NM/nn/de` identical to upstream (AS differs by design:
  expected-score N penalty, not the flat `--score-N`).

Real-data regression (shipped genomes, ~16 kb, real ACGT divergence):

```
for m in "-c --cs=long --MD" "-a --cs --MD" "-c -x asm20"; do
  diff <(./minimap2.upstream $m test/MT-human.fa test/MT-orang.fa | grep -v '^@PG') \
       <(./minimap2          $m test/MT-human.fa test/MT-orang.fa | grep -v '^@PG')
done   # -> no diff in any mode
```

## IUPAC-aware seeding — real-data result

On `Kmer2LTR_run.consensus.fa` (12,392 LTR consensus) vs
`gen5400000_final.fasta`, `-x asm20`:

| | stock | fork (aln-only) | fork (+IUPAC seeding) |
|---|---|---|---|
| consensus aligned | 7,676 | 7,458 | **7,803** |
| mean de (all primaries) | 0.0238 | — | **0.0004** |

- **+351 high-confidence elements vs stock** that stock + the aln-only
  fork both miss: median 393 bp, mapq 21, de 0.0027, **median 9.9% IUPAC
  content** (vs 1.9% overall) — i.e. exactly the IUPAC-dense copies the
  blind seeding couldn't anchor.
- 224 "lost vs stock" are the same junk churn (median 225 bp, mapq 0,
  de 0.04, 86 >5% div) — filtered out of any real soloLTR analysis.
- Discovery is capped below the optimistic projection sim (8,205) because
  the expected-score model honestly penalises degenerate columns (~ -1)
  rather than assuming a lucky concrete guess (+2). To push discovery
  further (at the cost of divergence-estimate conservatism) one could
  switch real-IUPAC columns to "compatible → ~full match" scoring — a
  one-function change in `mm_gen_iupac_mat`; not done (decision 1 stands).

## Notes / caveats

- `--score-N` now only governs the legacy ungapped short-read path /
  option validation; in the gapped path N is scored by the expected-score
  model (still a penalty). de/NM treatment of N is unchanged vs upstream.
- Minus-strand IUPAC works (`mm_idx_getseq_rev` + query revcomp use the
  symmetric `mm_comp_table`). When feeding a *query* that is itself an
  IUPAC string, complement it with an IUPAC-aware table (Python
  `maketrans("ACGT","TGCA")` does NOT complement R/Y/… — that is a
  caller-side concern, not minimap2).
- Out of scope (by decision): IUPAC-aware seeding; fractional identity;
  splice/jump logic beyond revcomp correctness.
- `minimap2.upstream` is a local diff artifact (git-ignored, not pushed).

## Resume this session

```
claude --resume "IUPAC-aware minimap2"
```
