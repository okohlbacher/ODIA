#!/usr/bin/env bash
# Run an adversarial review with the research vault wired in.
#
#   review_with_vault.sh -b BRIEF.md -t "topic words" [-s src1 -s src2 ...] [-o OUTDIR]
#
# Each reviewer gets the vault in the form it can actually consume:
#
#   kimi   -- absolute paths. It has a shell, but runs in a throwaway worktree where the
#             vault is absent (the vault lives outside the repo), so a relative path would
#             silently find nothing. Runs isolated: kimi has NO read-only mode and will
#             write to the tree it reviews.
#   codex  -- pasted note text plus any -s sources. codex 0.147.0 here cannot run shell
#             commands (its code-mode host binary was never installed); given a repo path
#             it returns an apology instead of a review, having read nothing.
#
# Identity is a PID, never a pgrep pattern: the waiting shell's own command line contains
# the pattern, so `pgrep -f "kimi -p"` matches itself and the wait never returns.

set -uo pipefail

export PATH="$HOME/.local/bin:$HOME/.kimi-code/bin:$PATH"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VAULT="${ODIA_VAULT:-/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/vault}"

BRIEF="" TOPIC="" OUT="" ; SRCS=()
while getopts "b:t:s:o:v:" o; do case $o in
  b) BRIEF=$OPTARG ;; t) TOPIC=$OPTARG ;; s) SRCS+=("$OPTARG") ;;
  o) OUT=$OPTARG ;; v) VAULT=$OPTARG ;;
  *) sed -n '2,12p' "$0" >&2; exit 2 ;;
esac; done

[[ -f "$BRIEF" ]] || { echo "need -b BRIEF.md" >&2; exit 2; }
[[ -d "$VAULT"  ]] || { echo "vault not found: $VAULT" >&2; exit 2; }
OUT="${OUT:-$(mktemp -d)}"; mkdir -p "$OUT"

# --- Phase 1: freeze the scope. Reviewers reconcile against what they are told. ---
SHA=$(git -C "$REPO" rev-parse --short HEAD)
DIRTY=$(git -C "$REPO" status --porcelain | wc -l)
if (( DIRTY )); then
  echo "WARNING: $DIRTY uncommitted path(s). Reviewers will be told, but committing is better." >&2
  git -C "$REPO" status --porcelain >&2
fi
echo "scope: HEAD $SHA, $DIRTY uncommitted" | tee "$OUT/scope.txt"

CTX_PATHS="$OUT/vault_paths.md"; CTX_PASTE="$OUT/vault_paste.md"
python3 "$REPO/scripts/vault_context.py" --vault "$VAULT" --topic "$TOPIC"           > "$CTX_PATHS"
python3 "$REPO/scripts/vault_context.py" --vault "$VAULT" --topic "$TOPIC" --paste   > "$CTX_PASTE"

hdr() { printf 'Review ONLY. Do not modify any file. HEAD is %s; %s uncommitted path(s).\n\n' "$SHA" "$DIRTY"; }

# --- kimi: paths, and an isolated worktree because it cannot be made read-only. ---
{ hdr; cat "$CTX_PATHS"; echo; cat "$BRIEF"; } > "$OUT/brief_kimi.md"
WT="$OUT/wt-kimi"
git -C "$REPO" worktree add -q --detach "$WT" HEAD || { echo "worktree failed" >&2; exit 1; }
( cd "$WT" && kimi -p "$(cat "$OUT/brief_kimi.md")" </dev/null ) >"$OUT/kimi.md" 2>&1 &
KPID=$!   # kimi 0.34.0 writes its ANSWER to stderr; 2>&1 is deliberate, not sloppy.

# --- codex: paste the vault and the sources; it has no shell. ---
{ hdr; cat "$CTX_PASTE"; echo
  if ((${#SRCS[@]})); then
    echo "########## PASTED SOURCES (line-numbered; cite these line numbers) ##########"
    for f in "${SRCS[@]}"; do echo "===== FILE: $f ====="; cat -n "$REPO/$f" 2>/dev/null || cat -n "$f"; done
  fi
  echo; cat "$BRIEF"
  echo; echo "You have NO working shell. Review from the pasted text only."
} > "$OUT/brief_codex.md"
codex exec -s read-only --skip-git-repo-check --color never \
      "$(cat "$OUT/brief_codex.md")" </dev/null 2>"$OUT/codex.log" >"$OUT/codex.md" &
CPID=$!

echo "kimi pid $KPID, codex pid $CPID -> $OUT"
wait $KPID; KRC=$?
wait $CPID; CRC=$?

# --- Phase 2 verification: prove kimi did not edit the tree it reviewed. ---
LEAK=$(git -C "$WT" status --porcelain | wc -l)
(( LEAK )) && { echo "!! kimi WROTE $LEAK path(s) in its worktree:" >&2; git -C "$WT" status --porcelain >&2; }
git -C "$REPO" worktree remove --force "$WT" 2>/dev/null

# codex fails closed and near-silently when its code-mode host is missing; catch it.
grep -qi "code-mode host\|could not read the required files" "$OUT/codex.md" 2>/dev/null &&
  echo "!! codex reviewed NOTHING (code-mode host missing) -- add more -s sources" >&2

printf 'kimi  rc=%s %s bytes -> %s\n' "$KRC" "$(wc -c <"$OUT/kimi.md")"  "$OUT/kimi.md"
printf 'codex rc=%s %s bytes -> %s\n' "$CRC" "$(wc -c <"$OUT/codex.md")" "$OUT/codex.md"
