# Design decisions

Running log. Each entry records what was decided and why, so the reasoning
survives the conversation it was made in. Open questions are listed at the end
and move up here as they are settled.

Evidence referenced as *(DIA-NN §x)* or *(OSW §x)* points at the two hand-off
documents in the project root, which are the design inputs for ODIA — see
`01-constraints.md` for what building the stack established.

---

## D1 — Phases are cut as thin end-to-end slices, not per pipeline stage

**Decided.** Each phase delivers something complete and usable; later phases
deepen it. The alternative was one phase per pipeline stage.

**Why:** the prior project's ±1% run-to-run ID variation (6,487 / 6,565 / 6,433)
was treated as a noise floor and used to size experiments. It turned out to be
two determinism bugs, neither in the classifier *(OSW §3.3)*, and several "no
measurable effect" verdicts sit inside that band. Determinism and phase
instrumentation are therefore preconditions, not deliverables — which requires
something running end to end early.

## D2 — Phase 1 is the assay library subsystem

**Decided.** Phase 1 covers the whole library path and can terminate there:

- read an existing library: TSV (DIA-NN-compatible) or Parquet
- generate from FASTA: digest → precursors → ONNX-predicted fragment
  intensities and iRT → transitions, with or without decoys
- write: TSV for interoperability, Parquet for speed
- stop-after-library, so the stage is independently usable

**Why:** it is the input every later phase needs, it is independently testable,
and it is where the prior project's largest memory term originated.

## D3 — ODIA owns its transition data structure

**Decided.** A lightweight, ODIA-specific representation rather than OpenMS's
`OpenSwath::LightTransition`.

**Why:** `LightTransition` carries `std::string transition_name`,
`std::string peptide_ref` and `std::vector<std::string> peptidoforms` per
transition row. At 78.6 M transitions that is the ~471 M allocations and the
+32.65 GB library load measured in the prior project *(OSW §2.2)*, of which ~88%
was allocator debris rather than live data. OpenMS has already applied this
reasoning once — its fragment-type enum comments that it "reduces memory from
~32 bytes (std::string) to 1 byte" — but has not carried it through the
remaining fields.

Proposed layout (struct-of-arrays, CSR: precursors sorted by m/z, each owning a
contiguous transition range):

```
precursors:   precursor_mz f64 | irt f32 | im f32 | charge u8 | decoy u8
              | modseq_id u32 | protein_id u32 | trans_begin u32 | trans_count u16
transitions:  product_mz | library_intensity f32
              | type u8 | ordinal u8 | charge i8 | loss u8 | flags u8
strings:      one arena + offset table; dictionary indices only
```

Two properties follow, and both matter beyond memory:

- sorting precursors by m/z makes each isolation window a **contiguous slice**
  rather than a gather, so the library can be paged instead of resident — the
  same shape as mzPeak's row-group pruning on the spectral side;
- Parquet input **keeps its dictionary**: read the dictionary once, store
  indices, never materialise a per-row `std::string`. This is the concrete
  reason the Parquet path is faster than TSV, rather than a vague claim.

**Measured, human proteome** (20,416 proteins; 3,987,909 precursors, 47,347,795
transitions): **759 MiB of arrays, 16.8 bytes per transition**, 1.20 GB peak RSS,
26 s. (Was 19.1 before decoys began sharing their target's sequence handle
rather than interning a mutated copy — see D7.)

Projected to the predecessor's library (78.6 M transitions) that is ~1.5 GB,
against **32.65 GB RSS / 7.26 GB live** there — roughly 22x on resident memory
and 4.8x on live data. The gap between those two ratios is the arena
fragmentation the allocation count caused, which is the part interning removes.

The figure was first reported as 16.5 bytes/transition. That was wrong:
`footprintBytes()` counted the character arena but not `StringArena::lookup_`,
one hash node per distinct string, nor unused block capacity — an understatement
of about 2.6x on a string-heavy library. Adversarial review caught it. The point
of D3 is that the memory claim is *measured*, so the measurement has to include
what is actually held.

## D4 — Use GPUs for ONNX inference when present, fall back to CPU

**Decided.** Execution provider selected at runtime: CUDA where available,
CPU otherwise. The same binary must work on both.

**Why:** predicting fragment intensities and iRT for a proteome-scale library is
the one embarrassingly parallel, arithmetically dense step in ODIA, and it is
the step most likely to be run on heterogeneous hardware.

Consequences to settle in Phase 1:

- No GPU exists on `ibminode05`, where this work is happening, so the CPU path
  is the one we can test here; the CUDA path needs a GPU node. The cluster
  almost certainly has one — `/ceph/ibmi/abi/opt/nvhpc` holds the NVIDIA HPC SDK
  with CUDA 12.8 and 13.1 — but Slurm cannot reach a controller from this node
  and the other nodes are not reachable over ssh, so this is unconfirmed.
- Provider registration must be **defensive**: attempt CUDA, catch, continue on
  CPU. A CUDA-enabled `onnxruntime` build that hard-fails when the driver is
  absent would make the binary unusable on exactly the node we develop on.
- conda-forge ships `onnxruntime-cpp` 1.26.0 in both CPU and CUDA builds; which
  one the environment carries is a packaging decision, not just a runtime one.

## D5 — `transition_name` is not stored

**Decided.** Synthesised on output from the peptide reference and the fragment
annotation; never held in memory.

**Why:** it is a derived identifier with ~78.5 M distinct values stored 78.5 M
times — the single largest string cost in the library *(OSW §2.2)*. The
annotation it is built from is itself reconstructible from
(type, ordinal, charge, loss), which we store in 4 bytes.

**Accepted risk:** round-trip fidelity is lost if an upstream producer encodes
something non-derivable in that field. The reader should therefore warn, once,
when an input's transition names do not match what ODIA would synthesise, so a
silent mismatch cannot go unnoticed.

## D6 — m/z is stored in 4 bytes as `uint32` fixed-point at 1e-5 Th

**Decided.** 4-byte m/z, using fixed-point rather than `float32`.

**Why:** both are 4 bytes, but they distribute error differently, and what
matters is *relative* precision because every tolerance in DIA is expressed in
ppm. Half-quantum error:

| m/z | `f32` | `u32` @1e-4 Th | `u32` @1e-5 Th |
|---:|---:|---:|---:|
| 150 | 0.051 ppm | 0.333 ppm | **0.033 ppm** |
| 500 | 0.031 ppm | 0.100 ppm | **0.010 ppm** |
| 2000 | 0.031 ppm | 0.025 ppm | **0.0025 ppm** |

`u32` at 1e-5 Th is at least as good as `f32` everywhere and an order of
magnitude better at the high end, while `u32` at 1e-4 Th is notably worse at low
m/z — fixed-point error is constant in absolute terms, so the scale has to be
chosen against the *bottom* of the range, not the top. It spans 0–42,949 Th,
far beyond any useful range.

Against a ~10 ppm extraction tolerance this contributes ~0.3% of the tolerance
width, and it stays well below the sub-ppm scale at which mass calibration fits
its residuals. Decode to `double` on use; the cost is one multiply.

`library_intensity` stays `f32` — it is relative and normalised, so relative
precision of 1e-7 is far more than the prediction warrants.

## D7 — Decoys by residue mutation, following DIA-NN

**Decided.** DIA-NN's approach: mutate a residue near each terminus via a fixed
substitution table *(DIA-NN §2.6)*. Fragment m/z are **recomputed from the
mutated sequence**, not shifted — shifting is only valid for ions spanning one
mutated residue, and the long b/y ions span both.

The decoy row stores the **target's** sequence, as DIA-NN does. The cost is that
a decoy's fragment m/z can no longer be checked against its own row — reproducing
them requires applying the substitution table, which thereby becomes an unwritten
part of the file format. Two reasons: the
substitution table is many-to-one, so distinct targets collide on a decoy
sequence and the synthesised identifiers stopped being unique; and since the
decoy inherits the target's precursor m/z by design, storing the mutated
sequence left the row internally inconsistent by up to 76 Th.

**Two things this obliges us to do**, because the DIA-NN document itself flags
this design as the weak point of its FDR model — the decoy keeps the target's
precursor m/z, its iRT and its library intensity pattern, so it is searched in
the same isolation window over the same RT range with the same expected spectrum
shape, and only the fragment masses differ:

1. **Keep the decoy generator pluggable.** Pseudo-reverse and shuffle must be
   selectable, so Phase 3 can measure this choice against entrapment rather than
   inherit it.
2. **Label symmetry is mandatory downstream.** Any filter or gate must be
   applied to targets and decoys by an identical criterion. If decoys are ever
   admitted via their target partner, every searched target has cleared a bar its
   decoy never faced, and since such criteria correlate with score the result is
   an *anti-conservative* FDR *(OSW §3.4)*.

## D8 — The library stage is a stop-point in the main tool

**Decided.** No separate library executable; the main tool can terminate after
the library stage and write it out.

**Why:** keeps one process and one data model, which is the point of D3 — every
process boundary would mean re-materialising the table we are working to keep
flat.

---

## Open

| # | Question | Blocks |
|---|---|---|
| O1 | Which ONNX models? None are on disk; `onnxruntime-cpp` 1.26.0 is packageable but the weights are not sourced. AlphaPeptDeep, Prosit, or something the prior OpenDIAlyzer carried? | Phase 1 generate path |
| O6 | No TSV/Parquet library exists on disk to test the reader against; only DIA-NN's binary `.speclib` (2.08 GB, in the bench directories). Export one with DIA-NN? | Phase 1 read path |
| O7 | Which TSV dialects must the reader accept? DIA-NN's own column names and the OpenSWATH set differ; DIA-NN normalises via synonym lists *(DIA-NN §2.1)*. Read both, write DIA-NN's? | Phase 1 read/write |
