#!/bin/bash
#refit training BestMods with sexFemale (male=0, female=1), sexEffect (female=-1, male=1) and
#pb(logAge_days) + pvc(logAge_days, by = sexEffect) coding,
#then compile a pdf comparing centile fans & sex-specific worm plots across all 4 codings
#run from outside code dir
#usage: ./code/submit_sex_coding.sh --list <file> [--pdf <output.pdf>] [--compile_only TRUE]
#  --list: one BestMod path per line (sexMale-coded training models, total-corrected or not),
#          optionally followed by a tab and its training csv (otherwise found from the BestMod's path)
#          e.g. find $(realpath cv_sample_A_train) -path "*_pbmods/model_objs/*_BestMod.rds" > bestmods.txt

compile_only="FALSE"

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --list) list="$2"; shift 2 ;;
    --pdf) pdf="$2"; shift 2 ;;
    --compile_only) compile_only="$2"; shift 2 ;;
    *) echo "Unknown parameter: $1"; exit 1 ;;
  esac
done

if [[ -z "$list" ]]; then
  echo "Missing --list"
  exit 1
fi

list=$(realpath $list)
pdf=${pdf:-${list%.*}_sex_coding.pdf}
pdf=$(realpath -m $pdf)
mkdir -p ./code/jobfiles/sex_coding

n=$(grep -v -c -e '^[[:space:]]*$' -e '^#' $list)
echo "$n models in $list"
echo "pdf: $pdf"

if [[ "$compile_only" == "TRUE" ]]; then
  sbatch ./code/subjob_compile_sex_coding.sh $list $pdf
  exit 0
fi

#refits (one array task per model), then compile once they've all finished (even if some failed)
array_id=$(sbatch --parsable --array=1-$n ./code/subjobs_sex_coding.sh $list)
echo "submitted refits: $array_id"
compile_id=$(sbatch --parsable --dependency=afterany:$array_id ./code/subjob_compile_sex_coding.sh $list $pdf)
echo "submitted pdf compile: $compile_id"
