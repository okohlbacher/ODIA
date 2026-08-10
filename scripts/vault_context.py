#!/usr/bin/env python3
"""Select vault notes relevant to a review topic and emit them as reviewer context.

The research vault lives OUTSIDE this repository (see ODIA_VAULT below) because it
carries competitor critique and is not distributable with the source.

Why this exists: on 2026-08-10 a design document went through two adversarial review
rounds carrying three factually wrong claims about DIA-NN's mass calibration. Neither
reviewer caught them, because both were handed the design and not the evidence -- while
`diann.cpp` sat unread in the reference directory. Reviewers get the vault now.

Two consumers, because the two CLIs differ in what they can reach:

  kimi   has a working shell, but runs in a throwaway git worktree where the vault is
         absent (it is not in the repo). It needs the ABSOLUTE PATH, which this prints.

  codex  0.147.0 on this box cannot run shell commands at all (its code-mode host binary
         was never installed), so it needs the note text PASTED into the prompt.

Usage:
    vault_context.py --topic "mass calibration ppm recalibration"      # ranked list
    vault_context.py --topic "..." --paste            # full text, for codex
    vault_context.py --topic "..." --paste --budget 60000 > ctx.md
    vault_context.py --always                         # just the must-read notes
"""

import argparse
import os
import re
import sys
from pathlib import Path

ODIA_VAULT = Path(os.environ.get(
    "ODIA_VAULT", "/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/vault"))

# Read regardless of topic. These are the notes that have repeatedly prevented a wrong
# conclusion rather than merely informed one; the provenance note is first because
# without it a reviewer will attribute the previous codebase's properties to this one.
ALWAYS = [
    "00-MOC/Vault provenance and the two-codebase caveat.md",
    "80-Measurement/Measurement traps in this project.md",
]

# Cheap stop list -- these words are in every note and carry no selectivity here.
STOP = set("""a an the and or of to in on for with without is are was were be been being
this that these those it its as at by from into than then there their they we our you your
not no any all both each more most other some such only own same so too very can will just
dia diann odia openswath ms1 ms2 note notes""".split())


def tokenize(text):
    return [w for w in re.findall(r"[a-z0-9_]+", text.lower())
            if len(w) > 2 and w not in STOP]


def score_note(path, terms):
    """Title matches count far more than body matches: a note *about* mass calibration
    beats one that mentions it in passing, and the vault's titles are descriptive."""
    try:
        body = path.read_text(errors="replace")
    except OSError:
        return 0.0, ""
    title = path.stem.lower()
    body_l = body.lower()
    score = 0.0
    for t in terms:
        if t in title:
            score += 10.0
        n = body_l.count(t)
        if n:
            # Diminishing returns, so one note repeating a word cannot dominate.
            score += min(n, 8) * 0.5
    return score, body


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--topic", default="", help="free text; words are matched against titles and bodies")
    ap.add_argument("--paste", action="store_true", help="emit full note text (for codex)")
    ap.add_argument("--budget", type=int, default=80000, help="max bytes of pasted text")
    ap.add_argument("--top", type=int, default=8, help="max topic-matched notes")
    ap.add_argument("--always", action="store_true", help="only the must-read notes")
    ap.add_argument("--vault", default=str(ODIA_VAULT))
    args = ap.parse_args()

    vault = Path(args.vault)
    if not vault.is_dir():
        sys.exit(f"vault not found: {vault}\nSet ODIA_VAULT or pass --vault.")

    chosen = []
    for rel in ALWAYS:
        p = vault / rel
        if p.is_file():
            chosen.append((999.0, p))
        else:
            print(f"# WARNING: must-read note missing: {rel}", file=sys.stderr)

    if args.topic and not args.always:
        terms = tokenize(args.topic)
        scored = []
        for p in sorted(vault.rglob("*.md")):
            if any(str(p).endswith(rel) for rel in ALWAYS):
                continue
            if "/.obsidian/" in str(p):
                continue
            s, _ = score_note(p, terms)
            if s > 0:
                scored.append((s, p))
        scored.sort(key=lambda x: (-x[0], str(x[1])))
        chosen.extend(scored[:args.top])

    print("# Research vault context")
    print()
    print(f"Vault root: `{vault}`")
    print()
    print("Conventions you must respect when using these notes:")
    print()
    print("- **Every claim is anchored** to a measurement or a citation. Where a note gives a")
    print("  number, that number was measured or cited -- prefer it over your own estimate, and")
    print("  say so if you think it is wrong.")
    print("- **Refuted claims are marked refuted rather than deleted.** A note saying REFUTED is")
    print("  telling you what was believed and why it was wrong. Do not re-propose it.")
    print("- **The \"ODIA\" in older notes is a PREVIOUS codebase**, not the repository under")
    print("  review. Read the provenance note first. Claims about DIA-NN, OpenSWATH, the")
    print("  literature and measurement discipline transfer; claims about \"ODIA\" mostly do not.")
    print()

    if not args.paste:
        print("Selected notes (read them yourself -- you have a shell):")
        print()
        for s, p in chosen:
            tag = "MUST READ" if s >= 999 else f"score {s:.1f}"
            print(f"- [{tag}] `{p}`")
        print()
        print(f"Everything else: `find '{vault}' -name '*.md'` (76 notes). The maps of content")
        print(f"are in `{vault}/00-MOC/`.")
        return

    total = 0
    print("The following notes are pasted in full because you have no shell.")
    print()
    for s, p in chosen:
        body = p.read_text(errors="replace")
        if total + len(body) > args.budget:
            print(f"<!-- omitted for budget: {p.relative_to(vault)} ({len(body)} B) -->")
            continue
        total += len(body)
        print(f"===== VAULT NOTE: {p.relative_to(vault)} =====")
        print(body.rstrip())
        print()
    print(f"<!-- {total} bytes of vault context -->", file=sys.stderr)


if __name__ == "__main__":
    main()
