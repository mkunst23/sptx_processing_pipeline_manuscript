# hash:sha256:5f86094afe77555ae3d49281062be4c9ce4249b845f4df9aad1584c8fe028a61
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/py-r:python3.10.12-R4.2.3-IRkernel-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG AWS_ACCESS_KEY_ID
ARG AWS_DEFAULT_REGION
ARG AWS_SECRET_ACCESS_KEY

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN pip install -U --no-cache-dir \
    matplotlib==3.9.2 \
    numpy==1.26.4 \
    pandas==2.2.2 \
    scanpy==1.10.2 \
    seaborn==0.13.2 \
    squidpy==1.6.0

ARG RCACHE_OPTIONS
RUN echo $RCACHE_OPTIONS >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site
RUN Rscript -e 'remotes::install_version("Rvcg", "0.23")' \
    && Rscript -e 'remotes::install_version("alphashape3d", "1.3.2")' \
    && Rscript -e 'remotes::install_version("anndata", "0.7.5.6")' \
    && Rscript -e 'remotes::install_version("ape", "5.8")' \
    && Rscript -e 'remotes::install_version("bezier", "1.1.2")' \
    && Rscript -e 'remotes::install_version("colorRamps", "2.3.4")' \
    && Rscript -e 'remotes::install_version("cowplot", "1.1.3")' \
    && Rscript -e 'remotes::install_version("data.table", "1.15.4")' \
    && Rscript -e 'remotes::install_version("doParallel", "1.0.17")' \
    && Rscript -e 'remotes::install_version("dplyr", "1.1.4")' \
    && Rscript -e 'remotes::install_version("ggplot2", "3.5.1")' \
    && Rscript -e 'remotes::install_version("googlesheets4", "1.1.1")' \
    && Rscript -e 'remotes::install_version("grid", "4.2.3")' \
    && Rscript -e 'remotes::install_version("gridExtra", "2.3")' \
    && Rscript -e 'remotes::install_version("gridGraphics", "0.5-1")' \
    && Rscript -e 'remotes::install_version("jpeg", "0.1-10")' \
    && Rscript -e 'remotes::install_version("klaR", "1.7-3")' \
    && Rscript -e 'remotes::install_version("randomcoloR", "1.1.0.1")' \
    && Rscript -e 'remotes::install_version("rearrr", "0.3.4")' \
    && Rscript -e 'remotes::install_version("reticulate", "1.38.0")' \
    && Rscript -e 'remotes::install_version("s2", "1.1.7")' \
    && Rscript -e 'remotes::install_version("stringr", "1.5.1")' \
    && Rscript -e 'remotes::install_version("tibble", "3.2.1")' \
    && Rscript -e 'remotes::install_version("tidyr", "1.3.1")' \
    && Rscript -e 'remotes::install_version("units", "0.8-5")'

RUN Rscript -e 'remotes::install_github( \
        "ttecon/ttr", \
        ref = "52354d2f051b7c31b2ad8d9cb24cccded370b220")'

COPY postInstall /
RUN /postInstall
