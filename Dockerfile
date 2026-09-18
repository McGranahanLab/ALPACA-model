# syntax=docker/dockerfile:1.7
#
# ALPACA — Allele-specific Phylogenetic Analysis of Copy-number Aberrations
# https://github.com/McGranahanLab/ALPACA-model
#
# Solvers included:
#   - SCIP  (default, open-source, tested by authors)
#   - GLPK  (open-source, tested by authors)
#   - Gurobi (via gurobipy; requires a licence at RUN time — see README.md)


ARG MICROMAMBA_VERSION=1.5.10

FROM mambaorg/micromamba:${MICROMAMBA_VERSION} AS builder

USER root

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        ca-certificates \
        git \
 && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
WORKDIR /src

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml


RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba install -y -n base -c conda-forge \
        r-base \
        r-jsonlite \
        r-optparse \
 && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER . /src/
RUN pip install --no-deps --no-cache-dir . \
 && python -c "import alpaca; import alpaca.solvers" \
 && alpaca --help >/dev/null


FROM mambaorg/micromamba:${MICROMAMBA_VERSION} AS runtime

LABEL org.opencontainers.image.title="ALPACA" \
      org.opencontainers.image.description="Allele-specific Phylogenetic Analysis of Copy-number Aberrations" \
      org.opencontainers.image.source="https://github.com/McGranahanLab/ALPACA-model" \
      org.opencontainers.image.documentation="https://github.com/McGranahanLab/ALPACA-model/blob/master/README.md" \
      org.opencontainers.image.licenses="MIT"

USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        ca-certificates \
        chromium \
 && rm -rf /var/lib/apt/lists/*
ENV KALEIDO_BROWSER_EXECUTABLE=/usr/bin/chromium

COPY --from=builder --chown=$MAMBA_USER:$MAMBA_USER /opt/conda /opt/conda

USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENV GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic

WORKDIR /work

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "alpaca"]
CMD ["--help"]
