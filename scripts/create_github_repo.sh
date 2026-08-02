#!/usr/bin/env bash
#
# Create the private GitHub repository and push this checkout to it.
#
# Requires an authenticated gh. This machine has no stored GitHub credentials,
# so authenticate first -- once, interactively:
#
#   source scripts/env.sh
#   gh auth login --hostname github.com --git-protocol https --web
#
# then run this script. It is safe to re-run: if the remote already exists it
# is reused rather than recreated.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${here}/.." && pwd)"
# shellcheck source=/dev/null
source "${here}/env.sh"

SLUG="${ODIA_GITHUB_SLUG:-okohlbacher/ODIA}"

cd "${repo_dir}"

if ! gh auth status >/dev/null 2>&1; then
  cat >&2 <<EOF
Not authenticated to GitHub. Run:

  source scripts/env.sh
  gh auth login --hostname github.com --git-protocol https --web

then re-run this script.
EOF
  exit 1
fi

if gh repo view "${SLUG}" >/dev/null 2>&1; then
  echo "==> ${SLUG} already exists; reusing it"
else
  echo "==> creating private repository ${SLUG}"
  gh repo create "${SLUG}" \
    --private \
    --description "OpenDIAlyzer: TOPP-compatible DIA extraction built on OpenMS 3.6 and mzPeak" \
    --disable-wiki
fi

if git remote get-url origin >/dev/null 2>&1; then
  git remote set-url origin "https://github.com/${SLUG}.git"
else
  git remote add origin "https://github.com/${SLUG}.git"
fi

echo "==> pushing"
git push -u origin main

echo "==> done: https://github.com/${SLUG}"
