# Dockerfile used to create a jupyter notebook for ImmunoPhenoDB loading
# Last updated 09/26/2023

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
ARG OWNER=jupyter
ARG BASE_CONTAINER=$OWNER/scipy-notebook:python-3.10.11
FROM $BASE_CONTAINER

LABEL maintainer="Jupyter Project <jupyter@googlegroups.com>"

# Fix: https://github.com/hadolint/hadolint/wiki/DL4006
# Fix: https://github.com/koalaman/shellcheck/wiki/SC3014
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

USER root

# R pre-requisites
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    fonts-dejavu \
    gfortran \
    gcc && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

USER ${NB_UID}

# R packages including IRKernel which gets installed globally.
# r-e1071: dependency of the caret R package
RUN mamba install --yes \
    'r-base' \
    'r-biocmanager' \
    'r-caret' \
    'r-crayon' \
    'r-devtools' \
    'r-e1071' \
    'r-forecast' \
    'r-hexbin' \
    'r-htmltools' \
    'r-htmlwidgets' \
    'r-irkernel' \
    'r-nycflights13' \
    'r-randomforest' \
    'r-rcurl' \
    'r-rmarkdown' \
    'r-rodbc' \
    'r-rsqlite' \
    'r-shiny' \
    'r-tidymodels' \
    'r-tidyverse' \
    'rpy2' \
    'unixodbc' && \
    mamba clean --all -f -y && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

RUN R -e 'BiocManager::install("SummarizedExperiment", ask=FALSE)' && \
    R -e 'BiocManager::install("SingleR", ask=FALSE)'

# Python Packages
RUN pip install \
    'mysql-connector-python==8.0.32' \
    'gtfparse==1.3.0' \
    'pyensembl==2.2.8' \
    'umap-learn==0.5.3' \
    'plotly==5.13.1' \
    'statsmodels==0.13.5' \
    'seaborn==0.12.2' \
    'pynndescent==0.5.10' \
    'numba==0.58.0' \
    'scanpy==1.9.5'
    
WORKDIR "${HOME}"