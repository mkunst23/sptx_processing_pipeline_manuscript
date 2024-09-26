# hash:sha256:8f5c769d53906640304d85b4a26b8eba9187dd476cea345ec66648dc72171b6f
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/py-r:python3.10.12-R4.2.3-IRkernel-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG AWS_ACCESS_KEY_ID
ARG AWS_DEFAULT_REGION
ARG AWS_SECRET_ACCESS_KEY

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev=7.81.0-1ubuntu1.13 \
        libx11-dev=2:1.7.5-1ubuntu0.3 \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y \
        r-vegan==2.6_6.1 \
    && conda clean -ya

RUN pip install -U --no-cache-dir \
    matplotlib==3.9.2 \
    numpy==1.26.4 \
    pandas==2.2.2 \
    scanpy==1.10.2 \
    seaborn==0.13.2 \
    squidpy==1.6.0

ARG RCACHE_OPTIONS
RUN echo $RCACHE_OPTIONS >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site
RUN Rscript -e 'remotes::install_version("BiocManager", "1.30.25")' \
    && Rscript -e 'remotes::install_version("DescTools", "0.99.57")' \
    && Rscript -e 'remotes::install_version("Rvcg", "0.23")' \
    && Rscript -e 'remotes::install_version("alphashape3d", "1.3.2")' \
    && Rscript -e 'remotes::install_version("anndata", "0.7.5.6")' \
    && Rscript -e 'remotes::install_version("ape", "5.8")' \
    && Rscript -e 'remotes::install_version("aplot", "0.2.3")' \
    && Rscript -e 'remotes::install_version("bezier", "1.1.2")' \
    && Rscript -e 'remotes::install_version("colorRamps", "2.3.4")' \
    && Rscript -e 'remotes::install_version("concaveman")' \
    && Rscript -e 'remotes::install_version("cowplot", "1.1.3")' \
    && Rscript -e 'remotes::install_version("dabestr", "2023.9.12")' \
    && Rscript -e 'remotes::install_version("data.table", "1.16.0")' \
    && Rscript -e 'remotes::install_version("doParallel", "1.0.17")' \
    && Rscript -e 'remotes::install_version("dplyr", "1.1.4")' \
    && Rscript -e 'remotes::install_version("entropy", "1.3.1")' \
    && Rscript -e 'remotes::install_version("ggforce", "0.4.2")' \
    && Rscript -e 'remotes::install_version("ggplot2", "3.5.1")' \
    && Rscript -e 'remotes::install_version("ggraph", "2.2.1")' \
    && Rscript -e 'remotes::install_version("ggrastr", "1.0.2")' \
    && Rscript -e 'remotes::install_version("ggridges", "0.5.6")' \
    && Rscript -e 'remotes::install_version("googlesheets4", "1.1.1")' \
    && Rscript -e 'remotes::install_version("graphlayouts", "1.1.1")' \
    && Rscript -e 'remotes::install_version("gridExtra", "2.3")' \
    && Rscript -e 'remotes::install_version("gridGraphics", "0.5-1")' \
    && Rscript -e 'remotes::install_version("igraph", "2.0.3")' \
    && Rscript -e 'remotes::install_version("jpeg", "0.1-10")' \
    && Rscript -e 'remotes::install_version("klaR", "1.7-3")' \
    && Rscript -e 'remotes::install_version("lsa", "0.73.3")' \
    && Rscript -e 'remotes::install_version("paletteer", "1.6.0")' \
    && Rscript -e 'remotes::install_version("philentropy", "0.8.0")' \
    && Rscript -e 'remotes::install_version("randomcoloR", "1.1.0.1")' \
    && Rscript -e 'remotes::install_version("rearrr", "0.3.4")' \
    && Rscript -e 'remotes::install_version("reshape", "0.8.9")' \
    && Rscript -e 'remotes::install_version("reticulate", "1.38.0")' \
    && Rscript -e 'remotes::install_version("rjson", "0.2.9")' \
    && Rscript -e 'remotes::install_version("s2", "1.1.7")' \
    && Rscript -e 'remotes::install_version("scales", "1.3.0")' \
    && Rscript -e 'remotes::install_version("sjmisc", "2.8.10")' \
    && Rscript -e 'remotes::install_version("stringr", "1.5.1")' \
    && Rscript -e 'remotes::install_version("tibble", "3.2.1")' \
    && Rscript -e 'remotes::install_version("tidyr", "1.3.1")' \
    && Rscript -e 'remotes::install_version("units", "0.8-5")'

RUN Rscript -e 'remotes::install_github( \
        "AllenInstitute/scrattch.vis", \
        ref = "e1a93534a6e7b65c93ff45da4b754c2dcc1e055c")' \
    && Rscript -e 'remotes::install_github( \
        "ttecon/ttr", \
        ref = "52354d2f051b7c31b2ad8d9cb24cccded370b220")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "BiocVersion", \
        "rhdf5" \
    ))' # Original versions: 3.16.0 2.42.1

COPY postInstall /
RUN /postInstall
