# Synthetic FreeSurfer 7 outputs for testing `sharing_code/freesurfer_to_datadict.R`

Three fake subjects. **The numbers are made up** — plausible in scale and
internally consistent (e.g. `CortexVol` equals the summed regional `GrayVol`,
`BrainSegVol - BrainSegVolNotVent` equals the summed ventricles), but not real
data. Use them only to check that the conversion runs and lands in the right
shape.

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

## Run it

Subject-directory mode:

```bash
Rscript sharing_code/freesurfer_to_datadict.R \
  --subjects-dir sharing_code/testdata_fs7 \
  --demographics sharing_code/testdata_fs7/demographics.csv \
  --id-col participant --age-col age_years --age-units years --sex-col sex \
  --study TEST --site-col site \
  --fs-version FS7_T1 \
  --out /tmp/test_datadict.csv
```

Group-table mode (same expected output, minus floating-point noise):

```bash
Rscript sharing_code/freesurfer_to_datadict.R \
  --lh-area sharing_code/testdata_fs7/tables/lh.aparc.area.tsv \
  --rh-area sharing_code/testdata_fs7/tables/rh.aparc.area.tsv \
  --lh-volume sharing_code/testdata_fs7/tables/lh.aparc.volume.tsv \
  --rh-volume sharing_code/testdata_fs7/tables/rh.aparc.volume.tsv \
  --lh-thickness sharing_code/testdata_fs7/tables/lh.aparc.thickness.tsv \
  --rh-thickness sharing_code/testdata_fs7/tables/rh.aparc.thickness.tsv \
  --aseg sharing_code/testdata_fs7/tables/aseg.volume.tsv \
  --demographics sharing_code/testdata_fs7/demographics.csv \
  --id-col participant --age-col age_years --age-units years --sex-col sex \
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
