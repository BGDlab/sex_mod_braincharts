#!/bin/bash

# Script to search for failed R jobs via error patterns in any .err files in code/jobfiles and subdirectories

if [ $# -eq 0 ]; then
    >&2 echo "provide one or more search strings for error file names, e.g. 'name_jobnumber' or multiple patterns separated by space"
    exit 1
fi

# Find all .err files in code/jobfiles and subdirectories matching any of the patterns
err_files=()
for pattern in "$@"; do
    while IFS= read -r file; do
        err_files+=("$file")
    done < <(find ./code/jobfiles -type f -name "*${pattern}*.err")
done

# Remove duplicates
err_files=( $(printf "%s\n" "${err_files[@]}" | sort -u) )

if [ ${#err_files[@]} -eq 0 ]; then
    echo "No matching .err files found for the given patterns."
    exit 0
    else
    echo "${#err_files[@]} files found, checking..."
fi

# Search for common error patterns once to avoid duplicate printing when multiple patterns hit the same lines
grep -RH --after-context=2 -E "Error in|Killed|halted|error" "${err_files[@]}" 2>/dev/null | awk '
BEGIN { filename = "" }
/^--$/ { next }
{
  if ($1 ~ /\.err:/) {
    newfile = substr($1, 1, index($1, ":")-1)
    if (newfile != filename) {
      filename = newfile
      print "\nFile: " filename
    }
    $1=""
    sub(/^:/,"")
  }
  line = $0
  # Highlight the error pattern line
  if (line ~ /Error in|Killed|halted|error/) {
    print "ERROR FOUND: " line
  } else {
    print line
  }
}'




