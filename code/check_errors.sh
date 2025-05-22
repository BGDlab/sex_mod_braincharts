#!/bin/bash

#script to search for failed R jobs via "Error in"

#check arg
if [ $# -eq 0 ]; then
    >&2 echo "provide search string for error file name, e.g. 'name_jobnumber'"
    exit 1
fi

#function

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

grep -RH --after-context=2 "Killed" ./code/jobfiles/*$1*.err  | awk '
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


grep -RH --after-context=2 "halted" ./code/jobfiles/*$1*.err  | awk '
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




