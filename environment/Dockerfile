# hash:sha256:c662a0a7442c00a831494527e5c117aa6fba7deb2488dc7b1c01e4371c956b7e
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
RUN Rscript -e 'remotes::install_version("anndata", "0.7.5.6")' \
    && Rscript -e 'remotes::install_version("cowplot", "1.1.3")' \
    && Rscript -e 'remotes::install_version("data.table", "1.15.4")' \
    && Rscript -e 'remotes::install_version("dplyr", "1.1.4")' \
    && Rscript -e 'remotes::install_version("ggplot2", "3.5.1")' \
    && Rscript -e 'remotes::install_version("googlesheets4", "1.1.1")' \
    && Rscript -e 'remotes::install_version("grid", "4.2.3")' \
    && Rscript -e 'remotes::install_version("gridExtra", "2.3")' \
    && Rscript -e 'remotes::install_version("randomcoloR", "1.1.0.1")' \
    && Rscript -e 'remotes::install_version("rearrr", "0.3.4")' \
    && Rscript -e 'remotes::install_version("reticulate", "1.38.0")' \
    && Rscript -e 'remotes::install_version("stringr", "1.5.1")' \
    && Rscript -e 'remotes::install_version("tibble", "3.2.1")' \
    && Rscript -e 'remotes::install_version("tidyr", "1.3.1")'

COPY postInstall /
RUN /postInstall
