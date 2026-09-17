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
solving the conda environment. The resulting image is roughly 4 GB (~2.5 GB of
conda env, plus Chromium for headless PDF plot rendering).

If your host uses systemd-resolved (Ubuntu 22.04+, Debian 12+, Fedora),
containers may not be able to resolve DNS on the default `bridge` network.
Pass `--network=host` to the build if you see `Temporary failure resolving`
errors:

```bash
docker build --network=host -t alpaca:latest .
```

## Run — open-source solvers

The image defaults to the Pyomo + SCIP backend, which needs no licence:

```bash
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$PWD:/work" \
    alpaca:latest run \
    --input_tumour_directory /work/examples/example_cohort/input/LTX0000-Tumour1 \
    --output_directory       /work/examples/example_cohort/output/LTX0000-Tumour1 \
    --plot_output_mode pdf \
    --genome_build hg19 \
    --solver pyomo --pyomo_solver scip
```

`--user "$(id -u):$(id -g)"` makes the container write output files as the
host user (the image otherwise runs as its bundled `mambauser`, uid 57439,
which cannot write into host-owned directories).

To run the full example script (which also exercises `input-conversion`,
`ancestor-delta`, `ccd`, `wgd`, and `plot-tumour`), override the entrypoint
to the micromamba activation shim so the conda env stays on `$PATH`:

```bash
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$PWD:/work" -w /work \
    --entrypoint /usr/local/bin/_entrypoint.sh \
    alpaca:latest bash examples/run_example.sh
```

The example script uses `--solver gurobi` (ALPACA's CLI default). The
example dataset is small enough to solve under the size-limited demo licence
bundled with `gurobipy`, so no licence file is required for the smoke test.
Pass `--solver pyomo --pyomo_solver scip` in your own runs if you don't have
a Gurobi licence and your data exceeds the demo size.

## Run — Gurobi backend

The `gurobipy` package is installed via `environment.yml`, but Gurobi itself
requires a licence for anything larger than its built-in demo size (~2000
variables). Provide your licence at run time — never bake it into the image.

### Floating (token server) licence

The most common academic and site-wide setup. The licence file is a small
text file containing `TOKENSERVER=<host>` (and optionally `PORT=<port>`,
default `41954`). Mount it into the container and make sure the container
can reach the token server:

```bash
docker run --rm \
    --user "$(id -u):$(id -g)" \
    --network host \
    -v "$PWD:/work" \
    -v "$HOME/gurobi.lic:/opt/gurobi/gurobi.lic:ro" \
    -e GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic \
    alpaca:latest run \
    --input_tumour_directory /work/examples/example_cohort/input/LTX0000-Tumour1 \
    --output_directory       /work/examples/example_cohort/output/LTX0000-Tumour1 \
    --solver gurobi
```

`--network host` is the simplest way to guarantee the container can reach
the token server on the campus/lab network. If you prefer to keep container
isolation, drop `--network host` and instead ensure the token server's host
and port are reachable from Docker's default bridge network (may require
firewall/DNS tweaks).

### Web License Service (WLS) / Named-User Academic

Works out of the box — no hardware fingerprint, no token-server network
requirement. Same command as above but without `--network host`.

### Node-locked (`grbgetkey`) licence

Node-locked licences fingerprint the licensed machine's MAC/hostname, which
the container does not share. Two workarounds:

1. Request a **WLS Academic** or **floating** licence from Gurobi (free
   swap for academics), then use the sections above.
2. Run the container with `--network host --hostname $(hostname)` and bind
   `/etc/hosts:/etc/hosts:ro` so the licence check sees the licensed
   machine. Fragile and non-portable — do not use on HPC.

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

Run with a Gurobi floating (token-server) licence — Singularity uses the
host network by default, so no extra flag is needed to reach the token
server:

```bash
singularity exec \
    --bind /path/to/gurobi.lic:/opt/gurobi/gurobi.lic:ro \
    --env GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic \
    alpaca.sif alpaca run --solver gurobi ...
```

WLS / Named-User Academic licences use the same command. Node-locked
licences generally do not work under Singularity for the same reason as
Docker — the container's MAC/hostname differ from the host — unless you
add `--hostname <licensed-host>`. Floating and WLS licences are the least
painful options on HPC.

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
- **`Permission denied` writing outputs** — pass `--user "$(id -u):$(id -g)"`
  so the container writes files as the host user.
- **`Kaleido requires Google Chrome to be installed`** — the image ships
  Chromium at `/usr/bin/chromium` and exports `KALEIDO_BROWSER_EXECUTABLE`.
  If you rebuild without Chromium, either install it in your derived image
  or run with `--plot_output_mode notebook` to skip the PDF renderer.
- **`Temporary failure resolving deb.debian.org` during build** — pass
  `--network=host` to `docker build` (see [Build](#build)).
