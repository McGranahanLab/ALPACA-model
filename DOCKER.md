# Running ALPACA in a container

This directory ships a `Dockerfile` that builds a self-contained ALPACA image
with all supported solvers (SCIP, GLPK, and — optionally — Gurobi) plus the R
runtime needed by `alpaca input-conversion`.

## Contents

- [Build](#build)
- [Run — open-source solvers](#run--open-source-solvers)
- [Run — Gurobi backend](#run--gurobi-backend)
- [Singularity / Apptainer](#singularity--apptainer)
- [Extending the image](#extending-the-image)
- [Troubleshooting](#troubleshooting)

## Build

From the repository root:

```bash
docker build -t alpaca:latest .
```

The build takes ~10–15 minutes on a typical laptop; most of the time is spent
solving the conda environment. The resulting image is roughly 2.5 GB.

## Run — open-source solvers

The image defaults to the Pyomo + SCIP backend, which needs no licence:

```bash
docker run --rm -v "$PWD:/work" alpaca:latest run \
    --input_tumour_directory /work/examples/example_cohort/input/LTX0000-Tumour1 \
    --output_directory       /work/examples/example_cohort/output/LTX0000-Tumour1 \
    --plot_output_mode pdf \
    --genome_build hg19 \
    --solver pyomo --pyomo_solver scip
```

To run the full example script (which also exercises `input-conversion`,
`ancestor-delta`, `ccd`, `wgd`, and `plot-tumour`), override the entrypoint:

```bash
docker run --rm -v "$PWD:/work" -w /work \
    --entrypoint bash alpaca:latest examples/run_example.sh
```

The example script uses `--solver gurobi` as ALPACA's CLI default; pass
`--solver pyomo --pyomo_solver scip` in your own runs if you don't have a
Gurobi licence available.

## Run — Gurobi backend

The `gurobipy` package is installed via `environment.yml`, but Gurobi itself
requires a licence for anything larger than its built-in demo size (~2000
variables). Mount your licence at run time — never bake it into the image.

```bash
docker run --rm \
    -v "$PWD:/work" \
    -v "$HOME/gurobi.lic:/opt/gurobi/gurobi.lic:ro" \
    -e GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic \
    alpaca:latest run \
    --input_tumour_directory /work/examples/example_cohort/input/LTX0000-Tumour1 \
    --output_directory       /work/examples/example_cohort/output/LTX0000-Tumour1 \
    --solver gurobi
```

### Licence-type caveats

| Licence type                     | Works in containers? | Notes                                                            |
| -------------------------------- | -------------------- | ---------------------------------------------------------------- |
| WLS / Named-User Academic (Web)  | Yes                  | No hardware fingerprint check; mount `gurobi.lic` as above.      |
| Floating (token server)          | Yes                  | Container must reach the token server (usually port `41954`).    |
| Node-locked (`grbgetkey`)        | Often no             | Fingerprints host MAC/hostname; the container's differ. See below. |

If you only have a node-locked academic licence, either:

1. Request a **WLS Academic** licence from Gurobi (free swap for academics), or
2. Run the container with `--network host --hostname $(hostname)` and mount
   `/etc/hosts`, so the licence check sees the licensed machine. This is
   fragile and not recommended for HPC.

## Singularity / Apptainer

Convert the Docker image to a SIF file. From a machine that has Docker:

```bash
# Push to a registry, then pull on the HPC:
docker tag alpaca:latest yourorg/alpaca:latest
docker push yourorg/alpaca:latest

singularity pull alpaca.sif docker://yourorg/alpaca:latest
```

Or convert locally from a running docker daemon:

```bash
singularity build alpaca.sif docker-daemon://alpaca:latest
```

Run with an open-source solver:

```bash
singularity exec alpaca.sif alpaca run \
    --input_tumour_directory /path/to/input \
    --output_directory       /path/to/output \
    --solver pyomo --pyomo_solver scip
```

Run with Gurobi, bind-mounting your licence:

```bash
singularity exec \
    --bind /path/to/gurobi.lic:/opt/gurobi/gurobi.lic:ro \
    --env GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic \
    alpaca.sif alpaca run --solver gurobi ...
```

Node-locked licences generally do not work under Singularity for the same
reason as Docker — the container's MAC/hostname differ from the host — unless
you use `--net --network=host` combined with matching `--hostname`. WLS
licences are the least painful option on HPC.

## Extending the image

### Add the `refphase` R package

The R conversion scripts (`extract_rephase_data.R`,
`convert_conipher_output.R`) only import `jsonlite` and `optparse` and read
their inputs with base R (`load`, `readRDS`). This works as long as the input
`.RData` / `.RDS` files contain plain lists/data-frames.

If your Refphase output stores S4 objects that require the `refphase` class
definitions to deserialise, extend the image:

```dockerfile
FROM alpaca:latest
USER root
RUN micromamba install -y -n base -c conda-forge r-remotes r-biocmanager \
 && R -e "remotes::install_bitbucket('schwarzlab/refphase', dependencies = TRUE)"
USER $MAMBA_USER
```

### Add CONIPHER

Same pattern; install from the McGranahanLab GitHub repository if your
CONIPHER `.tree.RDS` needs it to deserialise.

## Troubleshooting

- **`gurobi.GurobiError: No Gurobi license found`** — verify the licence file
  is mounted where `GRB_LICENSE_FILE` points, and that its permissions allow
  the container user (uid 1000 by default) to read it.
- **`Restricted license` warning in Gurobi output** — the demo licence is in
  use; either mount a real licence or switch to `--solver pyomo --pyomo_solver scip`.
- **SCIP failures on complex segments** — the authors observe SCIP fails on
  ~0.5% of segments (README §Solver selection). Fall back to Gurobi for those.
- **Very slow build** — most time goes into solving the conda environment;
  cache `~/.cache/pip` and use BuildKit (`DOCKER_BUILDKIT=1`).
