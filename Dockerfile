# syntax=docker/dockerfile:1.7
#
# ALPACA — Allele-specific Phylogenetic Analysis of Copy-number Aberrations
# https://github.com/McGranahanLab/ALPACA-model
#
# Multi-stage image that mirrors the conda install path from the README.
# The `runtime` stage produces the shipped image; the `builder` stage is
# discarded to keep the final size down.
#
# Solvers included:
#   - SCIP  (default, open-source, tested by authors)
#   - GLPK  (open-source, tested by authors)
#   - Gurobi (via gurobipy; requires a licence at RUN time — see DOCKER.md)
#
# The gurobipy wheel ships with a size-limited licence sufficient for the
# example dataset. Real workloads need a licence file bind-mounted at run
# time; see DOCKER.md for both Docker and Singularity instructions.
#
# R + jsonlite + optparse are included so `alpaca input-conversion` works
# on CONIPHER/Refphase/ASCAT RData/RDS inputs. The `refphase` R package
# itself is not installed by default (heavy Bioconductor deps); DOCKER.md
# shows how to extend the image if your RData contains refphase S4 objects.

ARG MICROMAMBA_VERSION=1.5.10

# ---------------------------------------------------------------------------
# Builder stage: solves the conda env and installs ALPACA from source.
# ---------------------------------------------------------------------------
FROM mambaorg/micromamba:${MICROMAMBA_VERSION} AS builder

USER root

# Ship a minimal build toolchain for pip's editable install and pyomo helpers.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        ca-certificates \
        git \
 && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
WORKDIR /src

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Install the environment from the repo's environment.yml, plus:
#   - kneed / matplotlib-base: listed in pyproject.toml but not environment.yml;
#     the recipe/meta.yaml runtime list has them, so we mirror that here.
#   - r-base + r-jsonlite + r-optparse: for the input-conversion R scripts
#     used by convert_conipher_output.R and extract_rephase_data.R.
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba install -y -n base -c conda-forge \
        kneed \
        matplotlib-base \
        r-base \
        r-jsonlite \
        r-optparse \
 && micromamba clean --all --yes

# Ensure subsequent RUN steps have the base env activated on PATH.
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install ALPACA itself from the source tree.
COPY --chown=$MAMBA_USER:$MAMBA_USER . /src/
RUN pip install --no-deps --no-cache-dir . \
 && python -c "import alpaca; import alpaca.solvers" \
 && alpaca --help >/dev/null

# ---------------------------------------------------------------------------
# Runtime stage: copy the solved env into a clean image with a fresh user.
# ---------------------------------------------------------------------------
FROM mambaorg/micromamba:${MICROMAMBA_VERSION} AS runtime

LABEL org.opencontainers.image.title="ALPACA" \
      org.opencontainers.image.description="Allele-specific Phylogenetic Analysis of Copy-number Aberrations" \
      org.opencontainers.image.source="https://github.com/McGranahanLab/ALPACA-model" \
      org.opencontainers.image.documentation="https://github.com/McGranahanLab/ALPACA-model/blob/master/DOCKER.md" \
      org.opencontainers.image.licenses="MIT"

USER root
# ca-certificates: needed for TLS out of the container.
# chromium: kaleido >=1.0 shells out to headless Chrome to render Plotly
#   figures for `--plot_output_mode pdf`. Without it, PDF plots fail and
#   the example script's final step falls back to a warning.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        ca-certificates \
        chromium \
 && rm -rf /var/lib/apt/lists/*
ENV KALEIDO_BROWSER_EXECUTABLE=/usr/bin/chromium

# Bring the fully solved conda env over from the builder.
COPY --from=builder --chown=$MAMBA_USER:$MAMBA_USER /opt/conda /opt/conda

USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Default location for user-supplied Gurobi licence (bind-mount at run time).
ENV GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic

# `/work` is the recommended bind-mount point for inputs/outputs.
WORKDIR /work

# `alpaca` is on PATH via the activated env; keep the entrypoint thin so users
# can also drop into a shell (`docker run --entrypoint bash ...`) if needed.
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "alpaca"]
CMD ["--help"]
