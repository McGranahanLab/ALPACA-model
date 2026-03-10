# alpaca_input_formatting

# Requirements

Required Python libraries:

pandas

rdata

# Running test

To run the script:

```bash
bash input_conversion.sh --tumour_id LTX1187 --copy_number_tool refphase --refphase_rData /path/to/refphase.RData --CONIPHER_tree_object /path/to/tree.RDS --output_dir /path/to/output
```

E.g.:

```bash
tumour_id="LTX1187"
refphase_rData="/nemo/project/proj-tracerx-lung/tracerx/_RELEASE/release_tx842/${tumour_id}/refphase/refphase_filt/${tumour_id}_refphase_run.RData"
CONIPHER_tree_object="/nemo/project/proj-tracerx-lung/tracerx/_RELEASE/release_tx842/${tumour_id}/mutation_trees/tumour_1/${tumour_id}_1.tree.RDS"
output_dir="/nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/tracerx800/input/test_cohort/${tumour_id}"


bash input_conversion.sh --tumour_id $tumour_id --copy_number_tool refphase --refphase_rData $refphase_rData --CONIPHER_tree_object $CONIPHER_tree_object --output_dir $output_dir
```

Optional:

- `--chromosome <value>` to process only one chromosome (e.g. `1`, `chr1`, `X`).

# Battenberg inventory format

When `--copy_number_tool battenberg` is used with `--battenberg_inventory`, the inventory file must contain one row per sample.

Required columns:

- `logr_segmented_path`: path to `<sample_name>.logRsegmented.txt.gz`
- `mutant_logr_path`: path to `<sample_name>_mutantLogR_gcCorrected.tab.gz`
- `heterozygous_baf_path` (or `baf_segmented_path`): path to `<sample_name>*BAFsegmented.txt.gz`
- `purity_ploidy_path`: path to `<sample_name>_battenbergA*purity_ploidy.txt*`

Optional columns:

- `sample` (or `sample_name`): sample ID. If missing, it is inferred from `logr_segmented_path`.
- `tumour_id` (or `tumor_id` / `case_id`): recommended for multi-tumour inventories.

Notes:

- Relative paths are resolved relative to the inventory file location.
- If multiple purity/ploidy files match, the one containing `default` in its file name is selected.
- BAF values are read from the  `BAFseg` column (if present); for headerless files, the 3rd column fallback is used.
- You can use a single inventory file for multiple tumours. If a `tumour_id` column exists, rows are filtered to the `--tumour_id` passed to `input_conversion.sh`.
- If there is no tumour column, all rows are read and non-matching samples are removed later by intersection with CONIPHER sample names. This works, but using a tumour column is safer and faster.
