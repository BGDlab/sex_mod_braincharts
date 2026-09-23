# Calculating Out-of-Sample Reference Scores

This folder holds the code for scoring new (out-of-sample) brain MRI data against
the reference models from Gardner et al. For each subject and phenotype it
returns a standardized score (z-score) and centile benchmarked on the brain chart models.
Any batch (site) effects of your data are estimated and removed along the way.

The pipeline has two steps:

1. **(Optional) `freesurfer_to_datadict.R`** reshapes FreeSurfer output into the
   correct CSV format (described in `data_dictionary.csv`).
2. **`oos_reference_scores.R`** scores that CSV against the reference models.

If your data is already in the format described in `data_dictionary.csv`, skip
straight to [Step 2](#step-2-score-your-data).

---

## Step 1 (optional): convert FreeSurfer output

`freesurfer_to_datadict.R` reads FreeSurfer stats and writes a CSV with one row per scan
that's scorable by `oos_reference_scores.R`, including the correct column names and units.

It can be used to compile and rename the MRI phenotypes (i.e. just FreeSurfer outputs), 
or, when passed a CSV containing demographics, will return all columns (i.e. covariates)
necessary for generating reference scores.

Data will ultimately be formatted to match the columns described in `data_dictionary.csv`, 
including:

- 68 Desikan-Killiany regions × {surface area `SA`, gray matter volume `GM`,
  cortical thickness `CT`}, e.g. `lh.SA.bankssts`
- 30 subcortical volumes, e.g. `SUBC.Left.Hippocampus`
- global measures: `GMV`, `sGMV`, `WMV`, `CSF`, `CBV`, `TBV`, `mean.CT`,
  `total.SA`
- covariates: `fs_version_SA`/`_CT`/`_GM`, `study_site`, `sexMale`,
  `logAge_days`

It only needs base R, with no extra packages.

> **WARNING:** this script assumes data is cross-sectional (one demographics row per
> participant), and does not include any QC filtering. Remove failed scans yourself before
> or after conversion.

### Input mode A: recon-all subject directories

Reads `<SUBJECTS_DIR>/<subj>/stats/{lh.aparc.stats, rh.aparc.stats, aseg.stats}`.

```bash
Rscript scoring_new_data/freesurfer_to_datadict.R \
  --subjects-dir /path/to/SUBJECTS_DIR \
  --demographics demo.csv --demo-id-col participant \
  --age-col age_years --age-units years --sex-col sex \
  --study MYSTUDY --site-col site \
  --fs-version FS7_T1 \
  --out my_datadict.csv
```

By default every folder under `--subjects-dir` that has a `stats/aseg.stats` is
converted. Pass `--subjects-list ids.txt` (one ID per line) to convert only some
subjects or to fix the row order. With a list, subjects that failed show up as
all-NA rows instead of being silently left out.

### Input mode B: `aparcstats2table` / `asegstats2table` group tables

```bash
Rscript scoring_new_data/freesurfer_to_datadict.R \
  --lh-area lh.area.tsv --rh-area rh.area.tsv \
  --lh-volume lh.vol.tsv --rh-volume rh.vol.tsv \
  --lh-thickness lh.thick.tsv --rh-thickness rh.thick.tsv \
  --aseg aseg.tsv \
  --demographics demo.csv --demo-id-col participant \
  --age-col age_years --age-units years --sex-col sex \
  --study MYSTUDY --site-col site \
  --fs-version FS6_T1 \
  --out my_datadict.csv
```

Make the tables with `aparcstats2table --hemi {lh,rh} --meas {area,volume,thickness}`
and `asegstats2table --meas volume`. Use either `--subjects-dir` or the tables,
not both.

### Demographics and labels

The demographics CSV is joined to the FreeSurfer subject IDs through
`--demo-id-col`, which defaults to `participant`. Without `--demographics`,
`sexMale`, `logAge_days` and `study_site` are left as NA and must be filled manually
prior to scoring.

| Option | What it does |
|---|---|
| `--age-col`, `--age-units` | Age at scan, in `days` (default), `weeks`, `months` or `years`. |
| `--age-type` | `birth` (default) adds 280 days of gestation to get post-conception age. Use `postconception` if the age column already is post-conception age. |
| `--ga-col`, `--ga-units` | Per-subject gestational age at birth, used instead of 280 days where it isn't NA. |
| `--sex-col` | Recoded to `sexMale` (1/0). `--male-values` / `--female-values` set which values count as male/female (defaults: `M,m,Male,male,MALE` and `F,f,Female,female,FEMALE`). Any other value becomes NA, with a warning. |
| `--study`, `--site-col` | `study_site` becomes `<study>_<site>`. Use a study label that won't clash with the training studies. |
| `--fs-version` | Fills `fs_version_SA`, `_CT` and `_GM`. Override them one at a time with `--fs-version-sa`/`-ct`/`-gm`. |

`--fs-version` should be one of the pipeline labels the models were trained on:
`FS7_T1`, `FS6_T1`, `FS7_T1T2`, `FS6_T1T2`, `FS5.3`, `FS6_T1FLAIR`,
`FS7_T1FLAIR`, `reconall-clinical`, `FSInfant`, `synthseg`, `synthseg-young`,
`synthseg-robust`, `Custom`. See `data_dictionary.csv` for details.

Details of how each FreeSurfer measure maps to a dictionary variable are in the
comment block at the top of `freesurfer_to_datadict.R`. Run it with `--help`
for the full list of options.

### Output

The CSV has a leading ID column (named by `--id-name`, default `participant`)
followed by the dictionary variables in dictionary order. Columns that are NA
for every subject are dropped. Columns that are missing for only some subjects
are reported once at the end of the run.

### Try it on the synthetic test data

`testdata_fs7/` contains fake outputs from FreeSurfer 7 for 3 synthetic subjects. 
See `testdata_fs7/README.md` for commands for both input modes. The
subject-directory version is:

```bash
Rscript scoring_new_data/freesurfer_to_datadict.R \
  --subjects-dir scoring_new_data/testdata_fs7 \
  --demographics scoring_new_data/testdata_fs7/demographics.csv \
  --demo-id-col participant --age-col age_years --age-units years --sex-col sex \
  --study TEST --site-col site \
  --fs-version FS7_T1 \
  --out scoring_new_data/testdata_fs7/test_output.csv
```

---

## Step 2: score your data

`oos_reference_scores.R` scores each phenotype in `--pheno_list` against the
reference models and adds z-score and centile columns to your data.

### Requirements

R with `dplyr`, `data.table`, `gamlss`, `devtools` and `remotes`. The script
also installs the development versions of [`gamlssTools`](https://github.com/BGDlab/gamlssTools) 
and [`gamlss2charts`](https://github.com/andy1764/gamlss2charts).

By default, helper functions and gamlss models are streamed from GitHub to save
memory and space on your local machine. To use the scripts offline, you must
save `oos_reference_scores_helper_funs.R` to the same directory as `oos_reference_scores.R`
and download the necessary gamlss models. Note that reference scores are calculated as 
an average across two models (each fit on a different split-half of the LBCC data, 
see Gardner et al for details); downloading only one set of models (i.e. only split A or B) 
may result in misleading reference scores. Pass the directory containing your local 
gamlss models to the script via `--model_dir`.

### Input data

`--df` must be a CSV whose columns follow `data_dictionary.csv`. Every row needs:

| Column | Notes |
|---|---|
| `study_site` | Batch variable (change it with `--batch`). New sites are fine: their effects are estimated from your reference rows and removed. |
| `sexMale` | Numeric 0/1 only. The script stops if it finds other values. |
| `logAge_days` | log10 of post-conception age in days. |
| at least one of `fs_version_SA`, `fs_version_CT`, `fs_version_GM` | Must use the training labels listed in Step 1. |
| the phenotype columns | Any phenotypes you want scored. |
| `TBV`, `mean.CT`, `total.SA` | Only needed with `--total TRUE` (see below). |

### Run it

```bash
Rscript scoring_new_data/oos_reference_scores.R \
  --df my_datadict.csv \
  --total FALSE \
  --ref_data "dx == 'CN'" \
  --out_file my_ref_scores.csv
```

Run with `--help` (or with no arguments) to print the usage.

### Options

Options can be written as `--key value` or `--key=value`. Hyphens and
underscores are interchangeable (`--min-ref` = `--min_ref`).

| Option | Required | Default | Description |
|---|---|---|---|
| `--df` | yes | | CSV of data to score (a local path or URL). |
| `--total` | yes | | `TRUE` scores against models that correct for total brain size; `FALSE` against uncorrected models. See below. |
| `--out_file` | yes | | Where to write the output CSV. Its folder must already exist. |
| `--ref_data` | | none | Which rows define the reference group for estimating batch effects. Defaults to all rows; give either a condition on columns of `--df`, e.g. `"dx == 'CN'"`, or a path to a CSV of reference rows. |
| `--min_ref` | | `75` | Minimum number of reference rows a site needs for its batch effect to be reliably estimated. Sites with fewer get NA scores. Default based on authors' testing. |
| `--batch` | | `study_site` | Column holding the batch (site) variable whose effects will be estimated and removed. |
| `--pheno_list` | | all 241 phenotypes (`all_phenos.txt` on GitHub) | To score just a subset, pass a text file of phenotypes to score, one phenotype column name per line. |
| `--model_dir` | | none (stream from GitHub) | Alternatively, can point to a folder containing brain chart models. |
| `--model_ref` | | `sharing` | Branch, tag or commit SHA of `BGDlab/sex_mod_braincharts` to stream models from. Ignored when `--model_dir` is set. Pin a tag or SHA for reproducible results. |

#### `--total`: correcting for brain size

With `--total TRUE`, each phenotype is modeled with a brain-size covariate:

- volumes (regional, subcortical and global) are corrected for `TBV`
- cortical thickness for `mean.CT`
- surface area for `total.SA`

Those columns must then be present and non-missing. Run the script twice, once
with `TRUE` and once with `FALSE`, to get both sets of scores.

#### `--ref_data`: estimating batch effects

When scoring out-of-sample data (i.e. subjects from study sites not seen in the
training data) offsets are estimated and removed before scoring. These estimates
can be done on the full group (i.e. all data in that batch, default), or, optionally, from 
a reference group (typically controls). Scores for all rows, patients included, are then 
relative to that reference group. Each site needs at least `--min_ref` reference
rows; otherwise scores are returned as NA. For more details, see [`gamlss2charts`](https://github.com/andy1764/gamlss2charts).

### Output

`--out_file` contains every column of the input `--df`, in the same row order,
plus two columns for each phenotype that was scored:

- `<pheno>_z`: standardized score (average of the split A and B models)
- `<pheno>_centile`: the corresponding centile, between 0 and 1

---
