#!/bin/bash
#submit one array job per config line (split x pheno), one array task per bootstrap
#with --rerun TRUE, only bootstraps without finished output (_uniform_sexdiffs.csv + _info.csv) are submitted
#can be run from any directory (paths below are absolute)
#usage: ./code/submit_boot_mods.sh --config <config file> [--n_boot 500] [--throttle 200] [--rerun TRUE/FALSE] [--lines "1 2"] [--dry_run TRUE]

n_boot=1000
throttle=200
rerun="FALSE"
lines=""
dry_run="FALSE"

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --config) config="$2"; shift 2 ;;
    --n_boot) n_boot="$2"; shift 2 ;;
    --throttle) throttle="$2"; shift 2 ;;
    --rerun) rerun="$2"; shift 2 ;;
    --lines) lines="$2"; shift 2 ;;
    --dry_run) dry_run="$2"; shift 2 ;;
    *) echo "Unknown parameter: $1"; exit 1 ;;
  esac
done

if [[ -z "$config" ]]; then
  echo "Missing --config"
  exit 1
fi

config=$(realpath $config)

BASE=/mnt/isilon/bgdlab_processing/Margaret/sex_mod_pvc
#log dir from subjobs_boot_mods.sh's #SBATCH --output/--error - slurm won't create it
mkdir -p $BASE/code/jobfiles/boot

#default: every line in the config
if [[ -z "$lines" ]]; then
  lines=$(awk -F'\t' '{print $1+0}' $config)
fi

#collapse a sorted list of ids into slurm ranges, e.g. "1 2 3 7" -> "1-3,7"
to_ranges() {
  awk 'NR==1 {s=$1; p=$1; next}
       $1==p+1 {p=$1; next}
       {printf "%s%s,", s, (p>s ? "-" p : ""); s=$1; p=$1}
       END {if (NR) printf "%s%s", s, (p>s ? "-" p : "")}'
}

for line in $lines
do
  save_path=$(awk -F'\t' -v L=$line '$1+0==L {print $5}' $config)
  pheno=$(basename $save_path)
  split=$(echo "$save_path" | grep -oE "cv_sample_[AB]_boot" | sed -E 's/cv_sample_|_boot//g')

  if [[ "$rerun" == "TRUE" ]]; then
    missing=()
    for id in $(seq 1 $n_boot); do
      prefix=$save_path/boot_$(printf "%03d" $id)
      if ! [[ -f ${prefix}_uniform_sexdiffs.csv && -f ${prefix}_info.csv ]]; then
        missing+=($id)
      fi
    done
    if [ ${#missing[@]} -eq 0 ]; then
      echo "line $line ($split $pheno): all $n_boot done, skipping"
      continue
    fi
    array_spec=$(printf '%s\n' "${missing[@]}" | to_ranges)
    echo "line $line ($split $pheno): ${#missing[@]} missing"
  else
    array_spec="1-${n_boot}"
  fi

  cmd="sbatch --job-name=boot_${split}_${pheno} --array=${array_spec}%${throttle} $BASE/code/subjobs_boot_mods.sh $config $line"
  echo "$cmd"
  if [[ "$dry_run" != "TRUE" ]]; then
    $cmd
  fi
done
