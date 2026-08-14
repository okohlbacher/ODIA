# Standard preamble for every ODIA adversarial-review brief

Prepend this verbatim. It exists because round 5 was run without it and two of
its conclusions were amended within the hour by measurements already in the
vault. See `vault/70-Adversarial/Adversarial review round 5.md` §6.

---

Much of this project's research lives in an Obsidian vault at
`/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/vault/` — measurements, benchmark results,
DIA-NN and OpenSWATH method notes, prior adversarial review rounds, calibration
findings, measurement traps, and the OpenMS API reference under `95-OpenMS-API/`
(queryable with `vault/.tools/ask.py`). **Read what is relevant before
answering.** Every claim there is anchored to a measurement or a citation, and
refuted claims are marked rather than deleted — the refutations are often the
most useful notes.

**Read `00-MOC/Vault provenance and the two-codebase caveat.md` FIRST.** The
"ODIA" in unmarked notes is a PREVIOUS attempt, not the current repository.
Claims about DIA-NN, OpenSWATH, the literature and measurement discipline
transfer intact; claims about that code do not.

If a claim of yours contradicts a vault measurement, say so explicitly and name
the flaw in how the measurement was taken — do not silently re-derive.

---

Notes on delivery:

- **codex** reviews from an empty workspace, so paste the relevant notes into the
  brief rather than trusting it to open paths, and run it at
  `-c model_reasoning_effort="max"`.
- **kimi** writes its answer to stdout and its banner/resume line to stderr at
  0.36+; capture the streams separately.
- **vibe** has returned 0 bytes on long briefs; check the byte count before
  treating an empty answer as a real one.
