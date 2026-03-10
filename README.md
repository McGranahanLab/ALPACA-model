# ALPACA
ALPACA is a computational method to infer allele- and clone-specific copy-number profiles of tumour clones from multi-sample bulk DNA sequencing. Read our [publication in Nature](https://www.nature.com/articles/s41586-025-09398-w) to learn more.
```
     _____ __    _____ _____ _____ _____
    |  _  |  |  |  _  |  _  |     |  _  |
    |     |  |__|   __|     |   --|     |
    |__|__|_____|__|  |__|__|_____|__|__|
    /\⌒⌒⌒/\
    (⦿   ⦿)
    ( 'Y' )
     (   )
     (   )
     (   )
     (~ ~~~~~~~~~~)
     ( ~ ~~   ~~  )
     ( ~  ~ ~  ~  )
     (~  ~~~~~   ~)
     │ │      │ │
     │ │      │ │
```
Repository containing core ALPACA code

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [ALPACA](#alpaca)
   * [Getting started](#getting-started)
      + [Installation](#installation)
      + [Testing installation](#testing-installation)
   * [Tutorial](#tutorial)
      + [Required inputs](#required-inputs)
         - [1. Fractional copy-numbers for each sample and each genomic segment](#1-fractional-copy-numbers-for-each-sample-and-each-genomic-segment)
         - [2. Confidence intervals associated with each allele-specific fractional copy-number](#2-confidence-intervals-associated-with-each-allele-specific-fractional-copy-number)
         - [3. Clone proportions table](#3-clone-proportions-table)
         - [4. Phylogenetic tree](#4-phylogenetic-tree)
      + [Example input file structure](#example-input-file-structure)
   * [Running ALPACA](#running-alpaca)
      + [Generating ALPACA input from BAM files](#generating-alpaca-input-from-bam-files)
      + [Running ALPACA using CONIPHER and Refphase outputs](#running-alpaca-using-conipher-and-refphase-outputs)
      + [Running ALPACA](#running-alpaca-1)
      + [Nextflow wrapper](#[nextflow-wrapper])
      + [Available options](#available-options)
         - [Solver selection](#solver-selection)

<!-- TOC end -->



<!-- TOC --><a name="getting-started"></a>
## Getting started

<!-- TOC --><a name="installation"></a>
### Installation

Start by cloning this repository:

```bash
git clone https://github.com/McGranahanLab/ALPACA-model.git
cd ALPACA-model
```

ALPACA is implemented in python and requires Linux or macOS. 

To install all the required dependencies use 'alpaca_conda.yml':

```bash
conda env create --name alpaca --file environment.yml
```

Next, install ALPACA with pip:

```bash
conda run -n alpaca pip install dist/*.whl
```

ALPACA ships with a Gurobi backend by default—please obtain a free academic license at [Gurobi](https://www.gurobi.com/academia/academic-program-and-licenses) if you plan to use it. Alternatively, you can select the new Pyomo backend to run ALPACA with open-source MILP solvers such as CBC or GLPK (see [Solver selection](#solver-selection)).

<!-- TOC --><a name="testing-installation"></a>
### Testing installation

To ensure that ALPACA works correctly after installation, activate the environment:

```bash
conda activate alpaca
```

then run:

```bash
bash examples/run_example.sh
```

This command should create output these output files:

```bash
ALPACA-model/examples/example_cohort/output/LTX0000-Tumour1
├── cn_change_to_ancestor.csv
└── final_LTX0000-Tumour1.csv
```

<!-- TOC --><a name="tutorial"></a>
## Tutorial

<!-- TOC --><a name="required-inputs"></a>
### Required inputs

This section describes inputs required by ALPACA. If you are using CONIPHER and Refphase as input to ALPACA, these input will be generated automatically see section [Running ALPACA using CONIPHER and Refphase outputs](#running-alpaca-using-conipher-and-refphase-outputs) below. Input for each tumour should be stored in a separate directory, i.e. each of the input tables should only contain data obtained from a single tumour.

<!-- TOC --><a name="1-fractional-copy-numbers-for-each-sample-and-each-genomic-segment"></a>
#### 1. Fractional copy-numbers for each sample and each genomic segment

These should be stored in a data frame with the following columns:

|segment|sample|cpnA|cpnB|tumour_id|
|--------|--------|--------|--------|--------|
|1_6204266_6634901|U_LTX0000_SU_T1.R1|3.2|2.0|LTX0000-Tumour1|
|1_6204266_6634901|U_LTX0000_SU_T1.R2|3.3|2.3|LTX0000-Tumour1|
|1_6204266_6634901|U_LTX0000_SU_T1.R3|3.4|2.0|LTX0000-Tumour1|

The table above shows the input for one genomic segment located on chromosome 1, starting at the base 6204266 and ending at 6634901 (encoded in the segment name as `<chr>_<start>_<end>`). Column 'sample' contains sample names of the tumour: this example contains 3 different samples (R1, R2 and R3) obtained from a single tumour (`U_LTX0000_SU_T1`). The sample names are arbitrary, but must be coherent within the entire input (including other input files). Fractional, allele-specific copy-numbers are stored in columns `cpnA` and `cpnB` and lastly column `tumour_id` stores the identifier of the tumour.

The segments are stored in the `ALPACA_input_table.csv` file

IMPORTANT

Pay special attention to underscore `_` character - it is used by ALPACA during file parsing and its usage must conform to the example pattern shown above. Do not use it in your tumour identifier.

<!-- TOC --><a name="2-confidence-intervals-associated-with-each-allele-specific-fractional-copy-number"></a>
#### 2. Confidence intervals associated with each allele-specific fractional copy-number

This table (called `ci_table.csv`) is similar to the ALPACA_input_table but contains lower and upper confidence intervals for each genomic segment.

|segment|sample|lower_CI_A|upper_CI_A|lower_CI_B|upper_CI_B|tumour_id|ci_value|
|--------|--------|--------|--------|--------|--------|--------|--------|
|10_38599060_42906137|LTX0000_SU_T1-R1|3.218|4.196|2.200|3.085|LTX0000-Tumour1|0.5|
|10_38599060_42906137|LTX0000_SU_T1-R2|1.468|1.695|2.703|2.977|LTX0000-Tumour1|0.5|

<!-- TOC --><a name="3-clone-proportions-table"></a>
#### 3. Clone proportions table

Table containing cellular prevalence of each clone in each sample, saved as comma separated file under the name `cp_table.csv`. This can be derived from cancer cell fractions (CCF), for example by subtracting the CCF values of children clones from CCF values of their parents. CONIPHER contains `compute_subclone_proportions` function which can be adapted to output the clone proportions.
This table contains an index column specifying clone names (which must match the name of clones in phylogenetic tree - see below) and one column for each sample. Proportions should sum to 1 in each sample, but small deviations from 1 are tolerated.

|clone|U_LTX0000_SU_T1.R1|U_LTX0000_SU_T1.R2|U_LTX0000_SU_T1.R3|
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

<!-- TOC --><a name="4-phylogenetic-tree"></a>
#### 4. Phylogenetic tree

Provide the SNV tree either as a JSON list-of-paths (`tree_paths.json`, preferred) or as a Newick file (`tree_paths.nwk`). If `tree_paths.json` is missing, ALPACA will automatically read `tree_paths.nwk` and convert it to the internal list-of-paths format described below.

JSON format: a list of arrays where each sub-array represents the phylogenetic path from the trunk (most recent common ancestor) to a terminal clone (leaf). For example, consider a simple tree with MRCA and three subclones. Subclones A and B are direct descendants of MRCA, and clone C is the child of clone B:

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

The same tree in Newick format would look like:

```
(A,(C)B)MRCA;
```

A more complex tree, with name of clones consistent with names used in the `cp_table.csv` above would look like this:

```json
[["clone12", "clone13", "clone14", "clone8"], ["clone12", "clone13", "clone14", "clone15"], ["clone12", "clone13", "clone16", "clone18", "clone1"], ["clone12", "clone19", "clone20"], ["clone12", "clone19", "clone21"]]
```

<!-- TOC --><a name="example-input-file-structure"></a>
### Example input file structure

Overall, for each tumour we expect the following files:

```bash
LTX0000-Tumour1
├── ALPACA_input_table.csv
├── ci_table.csv
├── cp_table.csv
└── tree_paths.json  (or tree_paths.nwk)
```

<!-- TOC --><a name="running-alpaca"></a>
## Running ALPACA

<!-- TOC --><a name="generating-alpaca-input-from-bam-files"></a>
### Generating ALPACA input from BAM files

For a tutorial on running ALPACA from BAM files, please see the tutorial at: [ALPACA pipeline](https://github.com/McGranahanLab/ALPACA-pipeline)

<!-- TOC --><a name="running-alpaca-using-conipher-and-refphase-outputs"></a>
### Running ALPACA using CONIPHER and Refphase outputs

We recommend using [CONIPHER](https://github.com/McGranahanLab/CONIPHER/blob/main/README.md) and [Refphase](https://bitbucket.org/schwarzlab/refphase/src/master/) outputs to generate input to ALPACA.

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

This is how your input directory for a single tumour should look like:

```bash
LTX0000-Tumour1
├── LTX0000-Tumour1-refphase-results.RData
└── LTX0000-Tumour1.tree.RDS
```

Using these two files, run the input-conversion command which should be available in your system after installing ALPACA.

E.g.:

```bash
tumour_id="LTX0000-Tumour1"
refphase_rData="examples/example_cohort/input/${tumour_id}/${tumour_id}-refphase-results.RData"
CONIPHER_tree_object="examples/example_cohort/input/${tumour_id}/${tumour_id}.tree.RDS"
conversion_output_dir="examples/example_cohort/input/${tumour_id}"

alpaca input-conversion \
    --tumour_id $tumour_id \
    --copy_number_tool refphase \
    --refphase_rData $refphase_rData \
    --CONIPHER_tree_object $CONIPHER_tree_object \
    --output_dir $conversion_output_dir
    
```

Make sure that your 'tumour_id' is the same as 'CASE_ID' in CONIPHER output and as 'patient_tumour' in Refphase output.

Converting the input might take a while and you will see this output while the program runs:

```bash
Tumour ID: LTX0000-Tumour1
Running input_conversion - it may take a few minutes
Argument 4 (examples/example_cohort/input/LTX0000-Tumour1/LTX0000-Tumour1-refphase-results.RData): Exists
Argument 6 (examples/example_cohort/input/LTX0000-Tumour1/LTX0000-Tumour1.tree.RDS): Exists
Argument 8 (examples/example_cohort/input/LTX0000-Tumour1): Exists
/Users/pp/miniforge3/envs/main/lib/python3.13/site-packages/alpaca/scripts/submodules/alpaca_input_formatting/input_conversion.sh
```

In some cases, CONIPHER outputs multiple phylogenetic tree. If there is a specific tree you would like to use, add CONIPHER_tree_index option:

selected_CONIPHER_tree_index=13

```bash
alpaca input-conversion \
    --tumour_id $tumour_id \
    --copy_number_tool refphase \
    --refphase_rData $refphase_rData \
    --CONIPHER_tree_object $CONIPHER_tree_object \
    --output_dir $conversion_output_dir \
    --CONIPHER_tree_index $selected_CONIPHER_tree_index
```

Default value of this argument is 1.

The following options for refphase conversion are available:
 
You can pass these options to the `alpaca input-conversion` helper (they are forwarded to the internal conversion scripts). Options marked (required) must be supplied; others have sensible defaults.

```bash
--heterozygous_SNPs_threshold <int> # (optional, default=5) Minimum supporting het-SNPs per phased segment used by convert_refphase.py: segments with fewer SNPS will be discarded.
--chromosome <str>               # (optional) Process only a single chromosome (e.g. 1, chr1, X)
--ci_value <float>                 # (optional, default=0.5) Confidence level for copy-numbers
--n_bootstrap <int>                # (optional, default=100) Number of bootstrap iterations used to calculate confidence intervals
--recalculate_not_updated_cns <0|1> # (optional, default=0) If 1, forces recalculation of copy-numbers for segments flagged as not-updated by refphase
--recalculate_updated_cns <0|1>   # (optional, default=0) If 1, forces recalculation of copy-numbers for segments flagged as updated by refphase
--recalculate_reference_cns  <0|1>   # (optional, default=0) If 1, forces recalculation of copy-numbers for segments flagged as 'reference' by refphase
```

You can choose the copy-number source tool:

```bash
--copy_number_tool <refphase|battenberg> # (optional, default=refphase)
```

For Battenberg input, use either an inventory file (recommended) or directory discovery:

```bash
--copy_number_tool battenberg
--battenberg_inventory <path>  # recommended
# OR
--battenberg_input_dir <path>  # optional auto-discovery mode
```

Battenberg inventory should contain one row per sample.
Required columns:
`logr_segmented_path`, `mutant_logr_path`, `heterozygous_baf_path` (or `baf_segmented_path`), `purity_ploidy_path`.
Optional columns:
`sample` (or `sample_name`) and `tumour_id` (or `tumor_id` / `case_id`).

You can use a single inventory file for multiple tumours. If a tumour column exists, rows are filtered to the `--tumour_id` provided to `alpaca input-conversion`.

In discovery mode, files are matched by sample prefix and purity/ploidy files are searched with:
`<sample_name>_battenbergA*purity_ploidy.txt*`.
BAF files are matched as `<sample_name>*BAFsegmented.txt*`, and values are read from the exact `BAFsegmented` column.
If multiple purity/ploidy files match, the file containing `default` in its name is preferred.

Notes:
- `--recalculate_not_updated_cns`. Refphase is using multiple samples to improve phasing. To do so, it first partitions genome of all samples (from the same tumour) into "consistent segments". After this step, each sample will have the same number of segments defined by the same breakpoints. This means, that in some sample, two segments in sequence might have the same copy-number, and yet be considered separate. Next, for each of the consistent segments, refphase attempts to perform phasing and updates the copy-numbers accordingly. This means that for "copy-number updated" segments the final copy-number will differ from the original copy numbers of the parent segment. This argument controls the fate of the other segmnets, i.e. the ones where phasing was not performed. We can either keep the original copy number (derived from all the SNPs on the parent segment), or update it using only the SNPs present on the consisent, non 'copy-number updated' segment.

- `--recalculate_updated_cns`. For segments where Refphase DID perform a copy-number update (i.e. the algorithm found evidence of allelic imbalance and reported a new copy-number), the reported copy-number may still differ slightly from what one would obtain by re-calculating using the same SNP-level equations and estimated purity/ploidy. When `--recalculate_updated_cns` is set to 1, the conversion step will re-calculate the copy-number for those "updated" segments from the underlying SNP-level data and then compute confidence intervals around that recalculated value. When set to 0 (the default), the conversion will keep the copy-number value provided by Refphase for updated segments and only compute confidence intervals centered on that reported value. Use this option when you want the conversion to re-derive copy-numbers from raw SNP evidence for segments that Refphase already updated.

- `recalculate_reference_cns` Recalculates the copy-number for segments marked as 'is_reference' True in Refphase. Default refphase behaviour is to recalculate and then round these copy numbers to nearest integers. Setting this option to '1' will trigger recalculation without the rounding, i.e. leaving the copy number for these segments in fractional state.

<!-- TOC --><a name="running-alpaca-1"></a>
### Running ALPACA

Once input is generated, ALPACA can be run with:

```bash
input_tumour_directory="examples/example_cohort/input/${tumour_id}"
output_directory="examples/example_cohort/output/${tumour_id}"
alpaca run \
   --input_tumour_directory "${input_tumour_directory}" \
   --output_directory "${output_directory}" \
   --genome_build hg19
```

Note: genome_build is only required for plotting, not for the copy-number inferenece.

Your input_tumour_directory should look like this:

```bash
LTX0000-Tumour1
├── ALPACA_input_table.csv
├── ci_table.csv
├── cp_table.csv
└── tree_paths.json  (or tree_paths.nwk)
```

but if you started from CONIPHER and Refphase, intermediary input files will also be present:

```bash
├── ALPACA_input_table.csv  <- ALPACA input
├── LTX0000-Tumour1-refphase-results.RData <- Refphase input
├── LTX0000-Tumour1.tree.RDS  <- CONIPHER input
├── ci_table.csv  <- ALPACA input
├── cp_table.csv  <- ALPACA input
├── phased_segs.tsv <- Intermediary input files
├── phased_snps.tsv <- Intermediary input files
├── purity_ploidy.tsv <- Intermediary input files
└── tree_paths.json  (or tree_paths.nwk) <- ALPACA input
```

Once ALPACA starts running, you will see the logo and the progress bar:

![ALPACA_run](resources/ALPACA_run.png "Image showing the output visible while running ALPACA")

ALPACA operates sequentially on each genomic segment. During this process, the `ALPACA_input_table.csv` will be decomposed into separate segment files for each segment so you will see additional `segments` directory in your input directory:

```bash
├── ALPACA_input_table.csv
├── LTX0000-Tumour1-refphase-results.RData
├── LTX0000-Tumour1.tree.RDS
├── ci_table.csv
├── cp_table.csv
├── phased_segs.tsv
├── phased_snps.tsv
├── purity_ploidy.tsv
├── segments
│   ├── ALPACA_input_table_LTX0000-Tumour1_10_38599060_42906137.csv
│   ├── ALPACA_input_table_LTX0000-Tumour1_10_42906138_45934831.csv
│   ├── ALPACA_input_table_LTX0000-Tumour1_10_45934832_74237001.csv
│   ├── ALPACA_input_table_LTX0000-Tumour1_10_74237002_135381927.csv
│   ├── ALPACA_input_table_LTX0000-Tumour1_10_95074_38406884.csv
│   ├── ALPACA_input_table_LTX0000-Tumour1_11_17317215_44596409.csv
└── tree_paths.json  (or tree_paths.nwk) <- ALPACA input
```

Solution for each segment is saved in the output directory. While the programme is running, you will see that output directory is populated with separate solution files:

```bash
ALPACA-model/examples/example_cohort/output/LTX0000
├── optimal_LTX0000-Tumour1_10_38599060_42906137.csv
├── optimal_LTX0000-Tumour1_10_42906138_45934831.csv
├── optimal_LTX0000-Tumour1_10_45934832_74237001.csv
├── optimal_LTX0000-Tumour1_10_74237002_135381927.csv
├── optimal_LTX0000-Tumour1_10_95074_38406884.csv
├── optimal_LTX0000-Tumour1_11_17317215_44596409.csv
├── optimal_LTX0000-Tumour1_11_193863_2164677.csv
├── optimal_LTX0000-Tumour1_11_2164678_3249658.csv
├── optimal_LTX0000-Tumour1_11_3249659_17317214.csv
└── optimal_LTX0000-Tumour1_11_44596410_49598207.csv
```

Once all the segments are done, these intermediary files are concatenated into single output files and deleted. The final output for a single tumour will contain two files:

```bash
ALPACA-model/examples/example_cohort/output/LTX0000-Tumour1
├── cn_change_to_ancestor.csv
└── ALPACA_output_LTX0000-Tumour1.csv
```

`ALPACA_output_LTX0000-Tumour1.csv` contains the clone-specific copy-numbers for both alleles:

|tumour_id|segment|clone|pred_CN_A|pred_CN_B|
|--------|--------|--------|--------|--------|
|LTX0000-Tumour1|2_41509_27282430|clone17|2|1|
|LTX0000-Tumour1|2_41509_27282430|clone20|2|1|
|LTX0000-Tumour1|2_41509_27282430|clone14|1|2|
|LTX0000-Tumour1|2_41509_27282430|clone9|1|1|
|LTX0000-Tumour1|2_41509_27282430|clone15|3|1|

`cn_change_to_ancestor.csv` contains the information on the copy number change between each clone and its ancestor:

|tumour_id|segment|clone|pred_CN_A|pred_CN_B|parent|parent_pred_cpnA|parent_pred_cpnB|cn_dist_to_parent_A|cn_dist_to_parent_B|
|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|
|LTX0000-Tumour1|2_41509_27282430|clone17|2|1|clone6|2|1|0|0|
|LTX0000-Tumour1|2_41509_27282430|clone20|2|1|clone4|3|1|-1|0|
|LTX0000-Tumour1|2_41509_27282430|clone14|1|2|clone1|1|2|0|0|
|LTX0000-Tumour1|2_41509_27282430|clone9|1|1|clone2|1|1|0|0|
|LTX0000-Tumour1|2_41509_27282430|clone15|3|1|clone11|2|1|1|0|

#### Plotting outputs and genome builds

Heatmap visualisations rely on chromosome length tables that match the genome build used to generate your tumour segments. ALPACA now downloads these tables from UCSC automatically and caches them under `~/.cache/alpaca/genomes`. Supply `--genome_build {hg19|hg38}` whenever plotting is enabled (the flag defaults to `hg19` for backwards compatibility, but you should set `--genome_build hg38` for GRCh38 inputs).

If you disabled plotting during `alpaca run` (for example with `--plot_output_mode none`), you can regenerate notebooks or PDFs later via:

```bash
alpaca plot-tumour \
   --input_directory "${input_tumour_directory}" \
   --output_directory "${output_directory}" \
   --plot_output_mode pdf \
   --genome_build hg38
```

The `plot-tumour` helper accepts the same `--chr-table` / `--mutation-table` overrides as the main run, so you can point it at a pre-downloaded chromosome table when working offline.

<!-- TOC --><a name="nextflow-wrapper"></a>

### Nextflow wrapper
Nextflow wrapper to facilitate parallel processing on HPC for large cohorts is available at:

https://github.com/McGranahanLab/ALPACA-nextflow


<!-- TOC --><a name="available-options"></a>

### Available options

```bash
--time_limit <int>
```

Time limit for single gurobi iteration, by defalt set to 60. For very complex cases, the iteration time might exceed 60 seconds, in which case the optimisation is stopped and current best solution is selected. Please not that ALPACA performs multiple iterations on each segment, this argument restricts time limit for single such iteration.

```bash
--cpu <int>

```

Number of available CPUs (default = 1). 

```bash
--overwrite_output <value>
```

```bash
--output_all_solutions <0|1>
```

When `--output_all_solutions` is enabled, ALPACA will save every computed solution for each processed segment into an `all_solutions/<segment>` subdirectory alongside the usual optimal output. For each segment this includes:

- a CSV with all model solutions (clone, pred_CN_A, pred_CN_B, complexity and an `elbow_offset` - how far from seleccted elbow each solution is),
- an elbow table CSV annotated with which complexity was selected by the knee-finding algorithm, and
- a simple PNG plot showing D_score vs allowed_complexity with the selected complexity marked.

Note: this option can produce a large number of files (one set per segment). We recommend using it only for selected segments when debugging or tuning model/optimization settings.

Controls whether ALPACA overwrites existing temporary files.
Allowed values: 0 (do not overwrite), 1 (overwrite, default).
In the default 'tumour' mode, ALPACA iterates sequentially over each segment, saving temporary .csv tables with solutions for each. It then concatenates all segment solutions into one final output file. On systems with time constraints, if ALPACA isn't allocated enough time, the run might be incomplete, resulting in only some segment solutions being present, but not the final file. In such situations, if the user restarts ALPACA, it will begin from scratch and overwrite all previously created files. To reuse these files, run ALPACA with the --overwrite_output 0 option. The default setting for this option is --overwrite_output 1 to prevent unintended reuse of temporary files.


```bash
--min_ci <float>
```

Minimum confidence-interval span.

If provided, ALPACA will ensure that all per-sample, per-allele confidence intervals (the columns `lower_CI_A`, `upper_CI_A`, `lower_CI_B`, `upper_CI_B` in `ci_table.csv`) have at least the specified span. Any interval that is tighter than `min_ci` will be expanded so that `upper_CI = lower_CI + min_ci` (while keeping `lower_CI >= 0`).

When `--min_ci` is used, ALPACA also writes a per-segment JSON report recording which samples/alleles had their CI expanded; the report files are written into the output directory under `segment_reports/` and are named `segment_report_<tumour_id>_<segment>.json` so downstream scripts or dashboards can consume the adjustments programmatically.

Example:

```bash
alpaca run --input_tumour_directory examples/example_cohort/input/LTX0000-Tumour1 \
           --output_directory examples/example_cohort/output/LTX0000-Tumour1 \
           --min_ci 0.05
```

This will enforce a minimum CI span of 0.05 for all allele-specific CI values and produce per-segment reports documenting any changes.

```bash
--extra_columns <col1> <col2> ...
```

List of additional columns to include in the output CSV files. By default, ALPACA outputs only the essential columns (`tumour_id`, `segment`, `clone`, `pred_CN_A`, `pred_CN_B`). Use this argument to request internal metrics or debugging information.

Supported columns:
- `complexity`: The tree complexity (number of events) of the solution.
- `gurobi_time`: The runtime (in seconds) of the Gurobi optimizer for the solution.
- `gurobi_gap`: The optimality gap of the solution (0.0 = optimal).
- `CI_score`: The confidence interval objective score.
- `D_score`: The distance objective score.

`gurobi_time`/`gurobi_gap` (and their objective-specific variants) are only populated when the Gurobi backend is active. When running via Pyomo they are returned as `-1` so that downstream tooling can detect the absence of solver-native metrics.

**Note on Multi-Objective Metrics:**
If you request `gurobi_time` or `gurobi_gap`, ALPACA will automatically include objective-specific metrics if they are available in the model run. For example, if running with both Distance (D) and Confidence Interval (CI) objectives (the default), requesting `gurobi_gap` will add:
- `gurobi_gap_D`: The gap specifically for the Distance objective.
- `gurobi_gap_CI`: The gap specifically for the CI objective.

Example:
```bash
alpaca run ... --extra_columns complexity gurobi_time gurobi_gap
```

```bash
--strict_gap <0|1>
```

Control the Gurobi MIP gap tolerance for improved reproducibility (default: 1, enabled).

When `--strict_gap 1` is set (default), ALPACA configures Gurobi to only stop optimization when:
1. The solver proves optimality (gap = 0), or
2. The time limit is reached

This is achieved by setting both `MIPGap = 0.0` and `MIPGapAbs = 0.0`, which forces the solver to find provably optimal solutions rather than stopping at "good enough" solutions within default tolerances.

**Why use strict_gap?**
- **Reproducibility**: Default gap tolerances can cause the solver to stop at different points across runs, leading to slightly different solutions
- **Determinism**: With strict_gap enabled, solutions are more deterministic (though time limits can still affect results)

**Important notes:**
- When strict_gap is enabled and segments hit the time limit before proving optimality, those segments will be reported in `run_gap_summary.csv`
- The `run_gap_summary.csv` file lists all segments where gap > 0, showing:
  - The final gap value
  - Why optimality wasn't reached (`time_limit`, `gap_tolerance`, or `other`)
  - Runtime and complexity information
- Set `--strict_gap 0` to use Gurobi's default gap tolerances (may finish faster but less reproducible)

Example:
```bash
alpaca run ... --strict_gap 1  # Enforce zero gap (default)
alpaca run ... --strict_gap 0  # Allow default gap tolerances
```

```bash
--solver <gurobi|pyomo>
```

Select the solver backend. `gurobi` (default) requires a licensed Gurobi installation and runs the model via Gurobi API. `pyomo` builds the ALPACA formulation with Pyomo and allows an open-source MILP solver via `--pyomo_solver`. Use this to run ALPACA in environments where Gurobi is not available. Note that advanced diagnostics such as `gurobi_time`/`gurobi_gap` columns are only produced when the Gurobi backend is used. We tested the software with `scip` and `glpk` solvers. Other solvers, such as `cbc` or `HiGHS` error for more complex cases.

```bash
--pyomo_solver <scip|glpk|gurobi|...>
```

When `--solver pyomo` is set, this flag selects the Pyomo solver plugin (default: `scip`). Ensure the corresponding command-line solver is installed and on your `$PATH`. Any solver supported by Pyomo for MILP/MIP (CBC, GLPK, CPLEX, Gurobi, etc.) can be used in principle, but only `scip` and `glpk` were tested. Features that rely on quadratic terms (e.g. positive variability penalties) currently require a solver with MIQP support; the Pyomo backend will raise an error if an unsupported configuration is detected.

```bash
--pyomo_solver_options key=value [key=value ...]
```

Optional list of additional parameters forwarded directly to the selected Pyomo solver. Use this to tune tolerances or enable solver-specific behaviour, e.g. `--pyomo_solver_options ratioGap=1e-5 maxNodes=1000`. Values are parsed as Python literals when possible (`1e-5`, `True`, `None`, etc.), otherwise passed as strings. Provided options override ALPACA's defaults.

```bash
--solver_logs <path>
```

Optional path (file or directory) where solver logs should be written regardless of backend. When unset, the Gurobi backend continues to honour `--gurobi_logs` while other backends remain silent unless their solver exposes its own logging options.

```bash
--plot_output_mode <notebook|pdf|none>
```

Control how ALPACA emits the final visualisations once segment optimisation finishes (default: `notebook`).

- `notebook`: writes `<tumour_id>_plots.ipynb` inside the tumour output directory. The notebook mirrors `dev/plotting/test_plotting.ipynb` and contains five cells (imports, config, heatmap A, heatmap B, CN changes) so you can explore interactive Plotly figures directly in VS Code or Jupyter.
- `pdf`: generates the static PDFs that older releases produced (`<tumour_id>_A_heatmap.pdf`, `<tumour_id>_B_heatmap.pdf`, `<tumour_id>_cn_changes_per_clone.pdf`).
- `none`: skips plotting entirely, which can be useful on headless servers that lack Plotly/Kaleido support or if you only need CSV outputs.

The option applies both to `alpaca run ...` (plots are produced automatically alongside `cn_change_to_ancestor.csv` and other summaries) and to `python -m alpaca.plotting ...` when you just want to re-render artefacts for an existing tumour directory.

<!-- TOC --><a name="solver-selection"></a>
#### Solver selection

ALPACA now supports multiple optimization backends:

- **Gurobi backend (default):** identical to previous releases, including lexicographic multi-objective support, IIS generation, and per-objective runtime/gap metrics. Requires a valid Gurobi + `gurobipy` installation. This is the recommended way of running the software!
- **Pyomo backend:** Allows running ALPACA with open-source solvers (SCIP, GLPK) or any other MILP solver exposed through Pyomo. This path is ideal for clusters without commercial licenses. Install the solver executable separately (for example, `conda install -c conda-forge glpk`) and then run `alpaca ... --solver pyomo --pyomo_solver glpk`.


Practical notes:


1. Strict gap enforcement maps to the closest option supported by the selected solver (e.g., `ratioGap=0` for CBC). Some solvers may ignore unknown options; ALPACA will still run but you should consult the solver logs if reproducibility is critical.
2. Selecting `scip` or `scipampl` automatically raises Pyomo's `ampl_command_timeout` to 10 seconds to prevent slow `scip --version` checks from failing availability detection. You can still override this (or any other) solver parameter through `--pyomo_solver_options`.
3. Objective-specific metrics (`gurobi_time_*`, `gurobi_gap_*`) and IIS reports are only available when the Gurobi backend is active. Pyomo runs populate the standard columns (`CI_score`, `D_score`, `complexity`, etc.) but leave Gurobi-specific columns at `-1` if requested via `--extra_columns`.
4. To compare solvers, run the same input twice, e.g. `alpaca run ... --solver gurobi` and `alpaca run ... --solver pyomo --pyomo_solver cbc`, and diff the resulting `optimal_*.csv` files.

#### Known issues:

1. `cbc` and `HiGHS` solvers fail in complex cases and in simple cases. Use `scip` or `glpk` instead.
2. Solver-native logs (via --solver_logs option) currently work only for Gurobi.
3. Most solver metrics work only for Gurobi.
4. On simulated dataset, SCIP and GLPK produce results only slightly worse to Gurobi (GLPK slightly better than SCIP), but both fail for a small number of segments (GLPK failed on 14 out of 90220 segment, SCIP failed on 490 out of 90220 segments). 
