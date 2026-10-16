# Synthetic FreeSurfer 7 outputs for testing `scoring_new_data/freesurfer_to_datadict.R`

Three fake subjects created by Claude Opus 4.0.

```
sub-00{1,2,3}/stats/lh.aparc.stats   34 DK regions, FS7 column layout
sub-00{1,2,3}/stats/rh.aparc.stats
sub-00{1,2,3}/stats/aseg.stats       45 FS7 aseg structures + # Measure header
demographics.csv                     participant, age_years, sex, site
tables/                              the same values as aparcstats2table /
                                     asegstats2table group tables
```

The aseg files use the FS7 names `Left-Thalamus` / `Right-Thalamus` (FS6 called
them `*-Thalamus-Proper`), which is one of the aliases the script resolves.

## Running `freesurfer_to_datadict.R`

Subject-directory mode:

```bash
Rscript scoring_new_data/freesurfer_to_datadict.R \
  --subjects-dir scoring_new_data/testdata_fs7 \
  --demographics scoring_new_data/testdata_fs7/demographics.csv \
  --demo-id-col participant --age-col age_years --age-units years --sex-col sex \
  --study TEST --site-col site \
  --fs-version FS7_T1 \
  --out scoring_new_data/testdata_fs7/test_output.csv
```

Group-table mode (same expected output, minus floating-point noise):

```bash
Rscript scoring_new_data/freesurfer_to_datadict.R \
  --lh-area scoring_new_data/testdata_fs7/tables/lh.aparc.area.tsv \
  --rh-area scoring_new_data/testdata_fs7/tables/rh.aparc.area.tsv \
  --lh-volume scoring_new_data/testdata_fs7/tables/lh.aparc.volume.tsv \
  --rh-volume scoring_new_data/testdata_fs7/tables/rh.aparc.volume.tsv \
  --lh-thickness scoring_new_data/testdata_fs7/tables/lh.aparc.thickness.tsv \
  --rh-thickness scoring_new_data/testdata_fs7/tables/rh.aparc.thickness.tsv \
  --aseg scoring_new_data/testdata_fs7/tables/aseg.volume.tsv \
  --demographics scoring_new_data/testdata_fs7/demographics.csv \
  --demo-id-col participant --age-col age_years --age-units years --sex-col sex \
  --study TEST --site-col site \
  --fs-version FS7_T1 \
  --out /tmp/test_datadict_tables.csv
```

Both should write 3 rows x 249 columns (the 248 data dictionary variables plus
the `participant` ID column).

## Things worth breaking on purpose

- delete a region row from an `lh.aparc.stats` file -> that region's SA/GM/CT go
  NA, and `mean.CT` / `total.SA` go NA for that subject
- rename `CerebralWhiteMatterVol` to `CorticalWhiteMatterVol` (FS 5.3 wording)
  -> `WMV` should still fill in
- drop `--demographics` -> `sexMale` and `logAge_days` come back NA with warnings
