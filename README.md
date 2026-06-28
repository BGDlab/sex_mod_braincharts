# Sex-Moderated Brain Charts

## Data Prep and QC
Nearly all data prep, filtering, etc, is in `lbcc_eda.Rmd`. This includes identifying and writing lists of
imaging-derived phenotypes (IDPs) to test, which are saved in `pheno_lists/`
Data from the Children's Hospital of Philadelphia was first retrieved using `build_your_cohort_dev_mg.Rmd`

## Sex Effect Significance Testing
Each analysis step is run using 3 scripts: 
- a 'config' script that writes a config file
- an R script
- a 'subjob' script that submits the R script as a job to the compute cluster using argument specified in the config file(s)

1. Prepare each split-half for model training: `code/config_cv_dfs.sh`, `code/prep_cv_dfs.R`, `code/subjobs_cv_dfs.sh`
```
#prep sex-moderated model dataframes
./code/config_cv_dfs.sh --log_age TRUE --total FALSE --logPheno FALSE
sbatch code/subjobs_cv_dfs.sh ./code/config_files/cv_sample_A_logPhenoFALSE_totalFALSE_logAgeTRUE_df_config.txt
sbatch code/subjobs_cv_dfs.sh ./code/config_files/cv_sample_B_logPhenoFALSE_totalFALSE_logAgeTRUE_df_config.txt

#prep dataframes with brain size covars for total-size-corrected models
./code/config_cv_dfs.sh --log_age TRUE --total TRUE --logPheno FALSE
sbatch code/subjobs_cv_dfs.sh ./code/config_files/cv_sample_A_logPhenoFALSE_totalTRUE_logAgeTRUE_df_config.txt
sbatch code/subjobs_cv_dfs.sh ./code/config_files/cv_sample_B_logPhenoFALSE_totalTRUE_logAgeTRUE_df_config.txt
```

2. Train sex-moderated models: `code/config_cv_mods.sh`, `code/fit_cv_mods.R`, `code/subjobs_cv_mods.sh`
```
./code/config_cv_mods.sh --total FALSE --log_age TRUE --sm "pb"
sbatch --array=1-478%200 /mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/subjobs_cv_mods.sh
```

3. Train sex-moderated models with correction for total brain size: `code/config_cv_mods.sh`, `code/fit_cv_total_mods.R`, `code/subjobs_cv_total_mods.sh`
```
./code/config_cv_mods.sh --total TRUE --log_age TRUE --sm "pb"
sbatch --array=1-478%200 /mnt/isilon/bgdlab_processing/Margaret/sex_mod_braincharts/code/subjobs_cv_total_mods.sh
```

4. Fit test models in alternate split-half: `code/config_cv_mods_test.sh`, `code/fit_cv_mods_test.R`, `code/subjobs_cv_mods_test.sh`
```
#sex-moderated models
#config
./code/config_cv_mods_test.sh --total FALSE --log_age TRUE
#submit
sbatch code/subjobs_cv_mods_test.sh ./code/config_files/cv_sample_A_totalFALSE_logAgeTRUE_test_config.txt
sbatch code/subjobs_cv_mods_test.sh ./code/config_files/cv_sample_B_totalFALSE_logAgeTRUE_test_config.txt
##############################################################
#sex-moderated, total-brain-size-corrected models
#config
./code/config_cv_mods_test.sh --total TRUE --log_age TRUE
#submit
sbatch code/subjobs_cv_mods_test.sh ./code/config_files/cv_sample_A_totalTRUE_logAgeTRUE_test_config.txt
sbatch code/subjobs_cv_mods_test.sh ./code/config_files/cv_sample_B_totalTRUE_logAgeTRUE_test_config.txt
```

5. Compile and visualize results, trajectory analyses: `review_test_models.Rmd`

### Sensitivity Analyses: Weighting by SurfaceHoles
1. Prepare each split-half for model training: `code/sensitivity_df_prep.Rmd`
2. Train models weighted by QC: `code/config_cv_mods_weighted.sh`, `code/fit_cv_mods_weighted.R`, `code/subjobs_cv_mods_weighted.sh`
```
#sex-moderated models
#config
./code/config_cv_mods_weighted.sh --total FALSE --log_age TRUE
#submit
sbatch code/subjobs_cv_mods_weighted.sh ./code/config_files/cv_sample_A_totalFALSE_logAgeTRUE_weighted_config.txt
sbatch code/subjobs_cv_mods_weighted.sh ./code/config_files/cv_sample_B_totalFALSE_logAgeTRUE_weighted_config.txt
##############################################################
#sex-moderated, total-brain-size-corrected models
#config
./code/config_cv_mods_weighted.sh --total TRUE --log_age TRUE
#submit
sbatch code/subjobs_cv_mods_weighted.sh ./code/config_files/cv_sample_A_totalTRUE_logAgeTRUE_weighted_config.txt
sbatch code/subjobs_cv_mods_weighted.sh ./code/config_files/cv_sample_B_totalTRUE_logAgeTRUE_weighted_config.txt
```

3. Fit test models in alternate split-half: `code/config_cv_mods_test_weighted.sh`, `code/fit_cv_mods_test_weighted.R`, `code/subjobs_cv_mods_test_weighted.sh`
```
#sex-moderated models
#config
./code/config_cv_mods_test_weighted.sh --total FALSE --log_age TRUE 
#submit
sbatch code/subjobs_cv_mods_test_weighted.sh ./code/config_files/cv_sample_A_totalFALSE_logAgeTRUE_weighted_test_config.txt
sbatch code/subjobs_cv_mods_test_weighted.sh ./code/config_files/cv_sample_B_totalFALSE_logAgeTRUE_weighted_test_config.txt
##############################################################
#sex-moderated, total-brain-size-corrected models
#config
./code/config_cv_mods_test_weighted.sh --total TRUE --log_age TRUE
#submit
sbatch code/subjobs_cv_mods_test_weighted.sh ./code/config_files/cv_sample_A_totalTRUE_logAgeTRUE_weighted_test_config.txt
sbatch code/subjobs_cv_mods_test_weighted.sh ./code/config_files/cv_sample_B_totalTRUE_logAgeTRUE_weighted_test_config.txt
```

4. Compile and visualize results: `review_centile_test.Rmd`

### Sensitivity Analyses: Restricting to Ages 2+
Refits the best cortical thickness models on participants age 2+ only, using the same training/test split structure.
1. Train age-2+ models in each split-half: `code/config_cv_mods_age2plus.sh`, `code/fit_cv_mods_age2plus.R`, `code/subjobs_cv_mods_age2plus.sh`
2. Fit test models in alternate split-half: `code/config_cv_mods_test_age2plus.sh`, `code/fit_cv_mods_test_age2plus.R`, `code/subjobs_cv_mods_test_age2plus.sh`
3. Compile and visualize results: `review_centile_test.Rmd`

## Case-Control Analyses of Normative Scores
1. Derive normative scores from each split-half's sex-moderated model: `code/subjobs_centile_calc.sh`, `code/centile_calc.R`, `code/config_centile_calc.sh`
```
#config
./code/config_centile_test.sh --total FALSE --log_age TRUE
#submit
sbatch code/subjobs_centile_test.sh ./code/config_files/cv_sample_A_totalFALSE_logAgeTRUE_centtest_config.txt
sbatch code/subjobs_centile_test.sh ./code/config_files/cv_sample_B_totalFALSE_logAgeTRUE_centtest_config.txt
```
2. Average cases' scores across split-halves and run analyses: `code/subjobs_centile_test.sh`, `code/centile_test.R`, `code/config_centile_test.sh`
```
#config
./code/config_centile_test.sh --total FALSE --log_age TRUE
#submit
sbatch --array=1-6 code/subjobs_centile_test.sh ./code/config_files/cv_mods_totalFALSE_logAgeTRUE_centext_config.txt
```

## Figures
Main manuscript figures were assembled using `format_figures.Rmd`. 
Nicely formatted centile fan plots were created using: `code/config_replot.sh`, `code/replot_centiles.R`, `code/subjobs_replot.sh`. These plots, together with the sex bias trajectories in median and variability for each phenotype, were compiled into a single supplemental PDF using `code/plot_all_phenotypes_combined.R` (train/test comparison pages followed by sex-difference trajectory pages).

Other scripts used for plotting, formatting, and viewing figures are `code/grab_pngs.R`, `code/plot_cv_brain.R`, and `code/subjob_grab_pngs.sh`.

## Misc
- `code/gamlss_fit_funs.R`: helper functions used to fit gamlss models
- `code/centile_cor.R`: check correlations between reference/normative scores derived from each split-half test model

### QC and Convergence Checks (`code/checks/`)
- `code/checks/check_errors.sh`: check errors from SLURM jobs
- `code/checks/config_check_convergence.sh`, `code/checks/check_convergence.R`, `code/checks/subjobs_check_convergence.sh`: stand-alone checks that models converged successfully
- `code/checks/check_correct_train_mod.sh`: audit `weighted_*` or `age2plus_*` train model objects against the original `BestMod.rds` selections to flag any phenotype mismatches; writes mismatched paths to `pheno_mismatches.txt`
- `code/checks/cleanup_mismatched_phenos.sh`: delete the mismatched train+test artifacts identified by `check_correct_train_mod.sh` so the affected phenotypes can be re-run with `--rerun TRUE` (dry-run by default; pass `--apply` to actually delete)
- `code/checks/test_rs_integral.R`: unit test for function calculating Riemann-Stieltjes integral `rs_integral()`
