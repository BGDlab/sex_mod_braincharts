#!/bin/bash
#source from subjob scripts after BASE is set: stops the job if the checkout at $BASE isn't on the
#branch this code expects, e.g. if someone switched branches there while jobs were queued.
#reads .git directly (no git call), so it works for any user and on nodes without git
#usage: source $BASE/code/check_branch.sh

EXPECTED_BRANCH=sex_pvc

gitdir=$BASE/.git
if [[ -f $gitdir ]]; then
  #worktree: .git is a file pointing to the real git dir
  gitdir=$(sed -n 's/^gitdir: //p' "$gitdir")
  [[ $gitdir != /* ]] && gitdir=$BASE/$gitdir
fi
branch=$(sed -n 's|^ref: refs/heads/||p' "$gitdir/HEAD" 2>/dev/null)

if [[ "$branch" != "$EXPECTED_BRANCH" ]]; then
  echo "ERROR: $BASE is on branch '${branch:-unknown}', expected '$EXPECTED_BRANCH' - not running" >&2
  exit 1
fi
echo "branch: $branch"
