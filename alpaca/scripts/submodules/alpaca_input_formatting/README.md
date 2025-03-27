# alpaca_input_formatting

# Requirements

Required R libraries:

rjson

argparse

CONIPHER

fst

# Running test

To run the script:

```bash
bash input_conversion.sh --tumour_id LTX1187 --refphase_rData /path/to/refphase.RData --CONIPHER_tree_object /path/to/tree.RDS --output_dir /path/to/output
```

E.g.:

```bash
tumour_id="LTX1187"
refphase_rData="/nemo/project/proj-tracerx-lung/tracerx/_RELEASE/release_tx842/${tumour_id}/refphase/refphase_filt/${tumour_id}_refphase_run.RData"
CONIPHER_tree_object="/nemo/project/proj-tracerx-lung/tracerx/_RELEASE/release_tx842/${tumour_id}/mutation_trees/tumour_1/${tumour_id}_1.tree.RDS"
output_dir="/nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/tracerx800/input/test_cohort/${tumour_id}"


bash input_conversion.sh --tumour_id $tumour_id --refphase_rData $refphase_rData --CONIPHER_tree_object $CONIPHER_tree_object --output_dir $output_dir
```
