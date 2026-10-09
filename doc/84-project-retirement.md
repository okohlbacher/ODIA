# 84 — Project retirement (2026-10-09)

ODIA development stopped on 2026-10-06 by the project owner's decision: on the public benchmark the best
engine reached about a quarter of DIA-NN's identifications at matched error, too far from the goal to continue
incrementally. This document records the final state so the work can be picked up or cited later.

## Final measured state (PXD047793 run 009, timsTOF Pro diaPASEF, 120 min)

Library: DIA-NN 2.0 prediction with MMTS-Cys (UniMod:39), 6,796,214 targets incl. entrapment (r = 0.151887),
13,429,845 precursors with decoys. Comparator: DIA-NN 2.0 on the same library, 93,900 precursors at q <= 0.01
(entrapment FDP 1.634%, 233 entrapments). ODIA numbers are N(233): targets accepted while tolerating the same
233 entrapments, i.e. at DIA-NN's own error level.

| arm | code | N(233) | share of DIA-NN | wall | peak RSS |
|---|---|---|---|---|---|
| reference | `r1/ms1par-pool` 0710221 | 21,916 | 23.3% | 16 h 13 min | 543 GiB |
| refit 1/K0 prior (fitted on run 010) | same binary, `-im_prior <spec>` | 24,186 | 25.8% | 15 h 58 min | 542 GiB |
| FRAGVEC sub-scores | `r1/integration` fb2562a, `-fragvec_scores final` | 22,881 | 24.4% | 15 h 56 min | 544 GiB |
| integration binary, all features off | `r1/integration` fb2562a | identical output to the reference | — | 16 h 01 min | 542 GiB |

DIA-NN: 40.5 min wall, 25.4 h CPU, 53.6 GB. ODIA's own q <= 0.01 counts (18-20k) sit at about 3% true FDP by
entrapment and should not be compared with DIA-NN's.

Findings of the last rounds worth keeping:

- **BlockPool retention** (`fix/blockpool-release`): the extractor reused memory blocks only at the exact
  same size and freed nothing until a pass ended; run 009's MS1 cycle drifts 1.593 -> 1.642 s, so sizes never
  recurred between chunks and every full run was OOM-killed at ~2.25 TB. Releasing the pool per chunk (plus
  best-fit reuse) bounds pass 1 at 534 GiB. Output-identical on the fixture and at full scale.
- **Frame-parallel MS1 build** (`explore/ms1-parallel`): the 15-19 h serial stage drops to ~1.5 h, output-identical.
- **Integration**: fragvec, Gate C hash calibration, scorer replay/export, per-frame TOF calibration, pool fix,
  all behind default-off flags, reproduce the reference byte-for-byte with every flag off (full scale).
- Where the gap is (from the project's earlier measurements): candidate generation is essentially complete;
  most missed DIA-NN identifications are peaks ODIA found but could not score past the FDR threshold.

## Code map (all on GitHub unless noted)

| ref | what |
|---|---|
| `main` | last reviewed mainline (v0.5.0 + chunked-extraction fix), this document, doc/83 step-back |
| `r1/integration` | best engine: ms1par + fragvec + gatec-hash + scorer-replay + tof-cal + pool fix; tag `r1-integration-final` |
| `r1/ms1par-pool` | the reference binary of the final numbers (ms1par + pool fix) |
| `fix/blockpool-release`, `fix/chunked-extraction` | the two memory/extraction fixes |
| `explore/*` | single-feature branches (each flag-gated, fixture-gated) |
| `backup/*`, `archive/*` | pre-rebase backups, the orphaned psink commit, and verbatim WIP snapshots |

Branches whose tree still carries in-house sample identifiers (pre-relabelling history) and WIP snapshots
containing them are not on GitHub; they are preserved, with everything else, in the git bundle kept with the
project's archive on the institute's Ceph storage (`ARCHIVE/ODIA-all-refs-20261009.bundle`).

## Records

Run registrations, logs, reads and analyses (text only; large intermediate data was deleted at retirement)
are archived with the project directory on the institute's storage (`shared/pxd/` for the PXD era:
R1B_ARMS.md, R1C_ARMS.md, MEM_PASS1.md, poolfix/POOLFIX.md, arms_done.log, arms/reads/).
