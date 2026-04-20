#!/bin/bash

# Script to search for failed R jobs via error patterns in any .err files in code/jobfiles and subdirectories

if [ $# -eq 0 ]; then
    >&2 echo "provide one or more search strings for error file names, e.g. 'name_jobnumber' or multiple patterns separated by space"
    exit 1
fi

# Build a single find call with all patterns OR'd together
find_args=()
for pattern in "$@"; do
    find_args+=(-o -name "*${pattern}*.err")
done
find_args=("${find_args[@]:1}")  # strip leading -o

# Collect unique files via associative array (no subshell sort needed)
declare -A seen
err_files=()
while IFS= read -r file; do
    if [[ -z "${seen[$file]}" ]]; then
        seen[$file]=1
        err_files+=("$file")
    fi
done < <(find ./code/jobfiles -type f \( "${find_args[@]}" \))

if [ ${#err_files[@]} -eq 0 ]; then
    echo "No matching .err files found for the given patterns."
    exit 0
else
    echo "${#err_files[@]} files found, checking..."
fi

grep -RH --after-context=2 -E "Error in|Killed|halted|error" "${err_files[@]}" 2>/dev/null | awk '
BEGIN { filename = "" }
/^--$/ { filename = ""; next }
{
  if (match($0, /^[^:]+\.err:/)) {
    newfile = substr($0, 1, RLENGTH - 1)
    $0 = substr($0, RLENGTH + 1)
    if (newfile != filename) {
      filename = newfile
      print "\nFile: " filename
    }
  }
  print $0
}'