#!/bin/bash

#script to search for failed R jobs via "Error in"
grep -RH --after-context=2 "Error in" ./code/jobfiles/*$1*.err  | awk '
BEGIN { filename = "" }
/^--$/ { next }
{
  if ($1 ~ /\.err:/) {
    newfile = substr($1, 1, index($1, ":")-1)
    if (newfile != filename) {
      filename = newfile
      print "\nFile: " filename
    }
    # remove filename from line
    $1=""
    sub(/^:/,"") #remove :
  }
  print $0
}'
