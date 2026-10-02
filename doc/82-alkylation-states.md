# 82 — Alkylation states: the audit, the defect, and the corrected baseline

2026-09-01. Trigger: the trace-forensics round found the IH1 DIA-NN reference
had searched cysteine UNMODIFIED. This document is the project-wide audit
(verdicts BY MASS, never by id strings: 400 sampled cysteine precursors per
library tested against unmodified vs +57.02146/Cys), the root cause, and the
required parameters going forward.

## The audit

| library | used by | verdict |
|---|---|---|
| dn_pred_cam.parquet | samelib arms + corrected DIA-NN ref | **CAM** (400/400, 0 free) |
| mix10k_acq_lib.parquet | mix10k fixture arms | **CAM** (400/400) |
| human_v2.parquet | ODIA own-lib e2e (IH1) | **CAM** (400/400) |
| astral_mix_lib.parquet | Astral arms (ODIA) | **CAM** (400/400) |
| corpus_lib.parquet | trace-model corpus era | **CAM** (400/400) |
| ih1_targets_osw.tsv | OpenSWATH IH1 | **CAM** (2,580 Carbamidomethyl annotations) |
| astral_osw.tsv | OpenSWATH Astral | **CAM-FREE — STALE**: 0 annotations in 130,692 rows; Mag-Net-era artifact (that sample had free thiols, doc/30); NOT comparable to the current neat sample |
| frozen DIA-NN IH1 reference (33,330) | the benchmark denominator | **CAM-FREE — THE DEFECT** (400/400 cysteine matches at unmodified mass) |

Every ODIA-side arm ever measured on IH1/mix/Astral searched the correct
chemistry. The defect is confined to the frozen IH1 DIA-NN reference (and the
stale Mag-Net-era OSW Astral library, already superseded by policy).

## Where it happened, and why

Phase: the REFERENCE'S LIBRARY GENERATION, before this project's second
attempt began. `shared/ref/diann/how_ih1_was_run.sh` searched
`full_lib.predicted.speclib`, a library predicted WITHOUT `--unimod4`.
DIA-NN's COMMAND LINE does not apply carbamidomethylation by default (the GUI
pre-checks it; the CLI does not), so the omission was silent. It was even
RECORDED — shared/ref/diann/README.md's own table says "searched CAM-free" —
but the number 33,330 was consumed everywhere as the reference regardless.
The Astral reference (re-cut 2026-08-17) used the correct recipe and is
unaffected.

Measured cost of the defect: the sample IS carbamidomethylated, so the
reference could not correctly match any of the ~24.8%-cysteine search space —
its accepted set is 1.9% cysteine (646, all at wrong chemistry) vs ODIA's
9.1%. All "recall of DIA-NN" numbers in docs 77-81 are measured against this
depressed, cysteine-blind denominator (their internal comparisons remain
valid — every arm shared the same reference — but the absolute recall levels
and the "36% missed" framing carry the caveat). 60.1% of the P1 port's
"complementary" IDs were simply outside the reference's search space.

## Required parameters (the record)

- **DIA-NN library generation**: `--fasta-search --predictor --gen-spec-lib
  --unimod4` (the Astral reference's recorded recipe). NEVER omit --unimod4
  on an alkylated sample; the CLI will not add it for you.
- **DIA-NN search**: the recorded flags (`--f <run.d> --lib <CAM lib>
  --fasta bench_reviewed_entrap.fasta --reannotate --threads 48
  --qvalue 0.01`); the .d folder, never the mzML (no MS2 spectra via mzML).
  DIA-NN 2.0 accepts the library as parquet.
- **ODIA / DIALibraryGenerator**: CAM is the shipped default — every audited
  ODIA library carries it. No action.
- **OpenSWATH**: IH1 conversion (ih1_targets_osw.tsv) is CAM-correct;
  astral_osw.tsv must be regenerated from a CAM source before any new Astral
  OSW number is quoted.
- **Samples**: IH1/IH2 (in-house diaPASEF) and Astral NEAT are
  carbamidomethylated preparations. The retired Mag-Net Astral file was the
  free-thiol exception and is gone from the benchmark (2026-08-17 policy).

## The corrected baseline

The corrected IH1 DIA-NN reference is being computed with the identical
recorded command and the CAM library (dn_pred_cam.parquet — the same library
ODIA's samelib arms search, making the corrected comparison same-library AND
same-chemistry). The corrected four-tool table lands in this document's
companion report when the run and the ODIA own-library P1 arm finish.
Until then: **the number 33,330 is SUPERSEDED for every purpose except
reproducing historical comparisons.**

## The corrected baseline (measured 2026-09-01, all arms CAM-correct)

**The corrected reference:** DIA-NN on IH1 with the CAM library = **37,334**
precursors at 1% FDR (vs the retired CAM-free 33,330; +12.0% = the cysteine
hole filled: 9.3% of accepted IDs are now cysteine-containing, vs 1.9%
before). Its entrapment FDP at the operating point is **1.13%** — DIA-NN is
essentially FDR-calibrated on this data. 5,539 genes / 4,943 protein groups.

**The table** (IH1 unless noted; matched entrapment FDP cells as t/e;
r = 0.1464 shared library, 0.1465 human_v2):

| arm | nominal q<=0.01 (true FDP at op) | FDP<=2% | FDP<=3% | FDP<=5% | recall vs 37,334 |
|---|---|---|---|---|---|
| **DIA-NN CAM (corrected ref)** | **37,334 (1.13%)** | ~37.3k | — | — | 100% |
| ODIA + P1 port, same library | 16,965 (2.21%) | 16,233 | 19,012 | 21,668 | 64.1% @ depth 37,334 |
| ODIA pre-port, same library | 13,735 (2.09%) | 14,894 | 17,152 | 20,626 | 63.9% |
| ODIA + P1 port, own library (human_v2) | 16,619 (2.63%) | 15,049 | 17,038 | 18,213 | — |
| ODIA pre-port, own library (2x2 arm) | 14,167 (**10.50%**) | 7,629 | 9,343 | 10,252 | — |
| OpenSWATH + pyProphet | 941 (frozen record) | — | — | — | — |
| *Astral neat (2nd instrument)* | DIA-NN 7,462 / ODIA 6,881 | flat vs control at every e>=10 cell | | | |

Reading notes, binding:
1. The samelib pair (rows 2-3) is the clean single-factor comparison; the
   own-lib pair is NOT (the pre-port 2x2 arm ran an older binary lineage with
   a 10.5% operating-point FDP — its matched cells show how much of the old
   own-lib number was miscalibration; the P1 own-lib arm nearly doubles it at
   matched FDP<=2%, 15,049 vs 7,629, but binary generations and the port are
   confounded in that pair).
2. ODIA's recall of the corrected reference is 64.1% at matched depth —
   essentially unchanged from the 64.0% against the old reference: ODIA
   recovers the newly-added cysteine peptides at the same rate as the rest
   (it always searched CAM), so the gap analysis of docs 77-81 carries over
   proportionally.
3. At a common true-FDP operating point (~1-2%): DIA-NN 37,334 vs ODIA
   16,233-16,965 — ODIA stands at ~43-45% of the corrected reference.
4. OSW's 941 is the frozen record (CAM-correct library); no rerun was
   performed. astral_osw remains stale (CAM-free, Mag-Net era).
