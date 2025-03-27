Repository with core ALPACA model code

TODO

create input with ci_tables - copy from simple pipeline
remove auth from the repo

Install package with:
pip install dist/*.whl


# Software requirements

# ALPACA
ALPACA is implemented in python - the easiest way to install all the required dependencies is to use 'alpaca_conda.yml':

```bash
conda env create --name alpaca --file alpaca_conda.yml
```
ALPACA is using Gurobi solver - please obtain free academic license before running the model at https://www.gurobi.com/academia/academic-program-and-licenses

# Full pipeline


# Required inputs

ALPACA requires the following multi-sample input for each patient:

## 1. Fractional copy-numbers for each sample and each genomic segment

These should be stored in a data frame with the following columns:

|segment|sample|cpnA|cpnB|tumour_id|
|--------|--------|--------|--------|--------|
|1_6204266_6634901|U_LTXSIM001_SU_T1.R1|3.2|2.0|LTXSIM001|
|1_6204266_6634901|U_LTXSIM001_SU_T1.R2|3.3|2.3|LTXSIM001|
|1_6204266_6634901|U_LTXSIM001_SU_T1.R3|3.4|2.0|LTXSIM001|

The table above shows the input for one genomic segment located on chromosome 1, starting at the base 6204266 and ending at 6634901 (encoded in the segment name as `<chr>_<start>_<end>`). Column 'sample' contains sample names of the tumour: this example contains 3 different samples (R1, R2 and R3) obtained from a single tumour (`U_LTXSIM001_SU_T1`). The sample names are arbitrary, but must be coherent within the entire input (including other input files). Fractional, allele-specific copy-numbers are stored in columns `cpnA` and `cpnB` and lastly column 'tumour_id' stores the identifier of the tumour.

The segments are stored in the `ALPACA_input_table.csv` file, but to facilitate parallel processing they are also stored separately in `segments` subdirectory, where each comma separated file contains data of only single segment. In this files segment identifiers are encoded in the file name and the file name must correspond to this pattern
`ALPACA_input_table_<tumour_id>_<segment_id>.csv`
for example:

```bash
ALPACA_input_table_LTXSIM001_11_65545995_65629325.csv
```

Pay special attention to underscore (`_`) character - it is used by ALPACA during file parsing and its usage must conform to the example pattern shown above.

## 2. Confidence intervals associated with each allele-specific fractional copy-number

This table (called `ci_table.csv`) is similar to the ALPACA_input_table but contains lower and upper confidence intervals for each genomic segment. If you are using Refphase to generate the input, running `run_conversion.sh` script will generate such table automatically (see Running ALPACA using CONIPHER and Refphase outputs below).

#TODO add example table

## 3. Clone proportions table

Table containing cellular prevalence of each clone in each sample, saved as comma separated file under the name `cp_table.csv`. This can be derived from cancer cell fractions (CCF), for example by subtracting the CCF values of children clones from CCF values of their parents. If you are using CONIPHER to generate the input, running `run_conversion.sh` script will generate such table automatically (see Running ALPACA using CONIPHER and Refphase outputs below). Alternatively, CONIPHER contains `compute_subclone_proportions` function which can be adapted to output the clone proportions.
This table contains an index column specifing clone names (which must match the name of clones in phylogenetic tree - see below) and one column for each sample. Proportions should sum to 1 in each sample, but small deviations from 1 are tolerated.
|clone|U_LTXSIM001_SU_T1.R1|U_LTXSIM001_SU_T1.R2|U_LTXSIM001_SU_T1.R3|
|--------|--------|--------|--------|
|clone1|0.0309|0.0006|0.1383|
|clone12|0.2810|0.0|0.0|
|clone13|0.0|0.0253|0.1112|
|clone14|0.1557|0.0|0.0021|
|clone15|0.0|0.0|0.1598|
|clone19|0.0|0.4785|0.2534|
|clone20|0.0202|0.4460|0.3176|
|clone21|0.0684|0.0|0.0174|
|clone8|0.4434|0.0495|0.0|
|clone16|0.0|0.0|0.0|
|clone18|0.0|0.0|0.0|

## 4. Phylogenetic tree

A json file (named `tree_paths.json`) containing the SNV tree structured encoded as an array of arrays - each of the sub-arrays represents the phylogenetic path from the trunk (most recent common ancestor) to a terminal clone (leaf). For example, consider a simple tree with MRCA and three subclones. Subclones A and B are direct descendants of MRCA, and clone C is the child of clone B:
```
       MRCA
       ├── A
       └── B
           └── C

```

Such tree would be encoded as following in ALPACA format:

```json
[['MRCA','A'],['MRCA','B','C]]
```
A more complex tree, with name of clones consistent with names used in the `cp_table.csv` above would look like this:

```json
[["clone12", "clone13", "clone14", "clone8"], ["clone12", "clone13", "clone14", "clone15"], ["clone12", "clone13", "clone16", "clone18", "clone1"], ["clone12", "clone19", "clone20"], ["clone12", "clone19", "clone21"]]
```

If you are using CONIPHER to generate the input, running `run_conversion.sh` script will generate such table automatically (see Running ALPACA using CONIPHER and Refphase outputs below).

## Example input file structure

Overall, for each tumour we expect the following files:

```bash
LTX000
├── ALPACA_input_table.csv
├── #TODO add CIs
├── cp_table.csv
├── segments
│   ├── ALPACA_input_table_LTXSIM001_11_65545995_65629325.csv
│   ├── ALPACA_input_table_LTXSIM001_16_81211548_81816818.csv
│   ├── ALPACA_input_table_LTXSIM001_1_6204266_6634901.csv
│   ├── ALPACA_input_table_LTXSIM001_22_46829252_47073142.csv
│   └── ALPACA_input_table_LTXSIM001_3_10077023_10191734.csv
└── tree_paths.json

```

# Running ALPACA from BAM files

# Running ALPACA using CONIPHER and Refphase outputs

We recommend using CONIPHER and Refphase outputs to generate input to ALPACA.
CONIPHER output directory should contain a 'tree object' for each patient. This object is save as RDS file with the following name: <CASE_ID>.tree.RDS

Refphase aggregates all outputs in one object:

```R
results <- refphase(refphase_input)
```

Make sure that this entire object is saved as a single RData object, for example by adding this code to your Refphase script:

```R
results <- refphase(refphase_input)
save(results, file = paste0(patient, "-refphase-results.RData"))
```

Using these two files, run the 'run_conversion.sh' script located in submodules/alpaca_input_formatting.

E.g.:

```bash
tumour_id="LTX000"
refphase_rData="examples/simple_pipeline/${tumour_id}/refphase/${tumour_id}-refphase-results.RData"
CONIPHER_tree_object="examples/simple_pipeline/${tumour_id}/conipher/output/${tumour_id}.tree.RDS"
output_dir="examples/simple_pipeline/${tumour_id}/alpaca/input"

input_conversion --tumour_id $tumour_id --refphase_rData $refphase_rData --CONIPHER_tree_object $CONIPHER_tree_object --output_dir $output_dir
```

Make sure that your 'tumour_id' is the same as 'CASE_ID' in CONIPHER output and as 'patient_tumour' in Refphase output.

Once input is generated, ALPACA can be run with:
