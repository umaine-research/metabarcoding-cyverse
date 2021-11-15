FROM cyversevice/rstudio-verse:latest

MAINTAINER Laura Jackson, laura.jackson@maine.edu"

ARG BUILD_DATE
ARG VCS_REF
ARG VERSION
LABEL org.label-schema.build-date=$BUILD_DATE \
      org.label-schema.name="RStudio Verse" \
      org.label-schema.description="Built from Rocker-Project RStudio Verse, additional depends for CyVerse K8s workbench" \
      org.label-schema.url="https://cyverse.org" \
      org.label-schema.vcs-ref=$VCS_REF \
      org.label-schema.vcs-url="e.g. https://github.com/cyverse-vice/rstudio-verse" \
      org.label-schema.vendor="CyVerse" \
      org.label-schema.version=$VERSION \
      org.label-schema.schema-version="1.0.0"

# Install Tools/Ubuntu Setup
#RUN apt-get update && apt-get install -y apt-utils

# Install R packages
USER root
RUN R -e 'install.packages(c("BiocManager", "zetadiv", "deblur", "decontam", "R.utils"), repos = "https://cloud.r-project.org/", dependencies = TRUE)'

# Set Biocmanager v3.11 and Install dependencies
RUN R -e 'BiocManager::install(version = "3.14", ask = FALSE)'
RUN R -e 'BiocManager::install("dada2", version = "3.14")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("DECIPHER")'
RUN R -e 'BiocManager::install("phyloseq")'
RUN R -e 'BiocManager::install("ShortRead")'


# RStudio Setup
USER rstudio
WORKDIR /home/rstudio

# Install Miniconda3
ENV PATH="/home/rstudio/miniconda3/bin:${PATH}"
ARG PATH="/home/rstudio/miniconda3/bin:${PATH}"
RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /home/rstudio/.conda
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

# Update Miniconda
RUN conda init bash \
  && . ~/.bashrc \
  && conda update -n base -c defaults conda \
  && conda config --add channels defaults \
  && conda config --add channels bioconda \
  && conda config --add channels conda-forge

# Install Cutadapt, FASTQC
RUN conda install -y cutadapt fastqc \
  && echo 'alias cutadapt = "conda run --prefix '/home/rstudio/miniconda3/envs/cutadaptenv' cutadapt"'

COPY profile .profile

USER root
ENTRYPOINT ["/usr/local/bin/run.sh"]
