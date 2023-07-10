FROM harbor.cyverse.org/vice/jupyter/datascience

USER root 

# Add repository for updated R version
RUN sudo apt-get update && apt-get install -y software-properties-common
RUN sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

# Import GPG key for R repository
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9

# Update and clean package manager
RUN apt-get update --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Update and install R
RUN sudo apt-get update && apt-get install -y r-base

# Set CRAN mirror to a specific location
RUN echo 'options(repos = "https://cloud.r-project.org")' >> /opt/conda/lib/R/etc/Rprofile.site

# Update and install packages
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'

RUN R -e 'BiocManager::install(version = "3.16")'
RUN R -e 'BiocManager::install(c("DESeq2", "bambu"))'

USER jovyan