FROM harbor.cyverse.org/vice/jupyter/rstudio:latest
RUN conda install -c conda-forge mamba -y
RUN mamba install --quiet --yes tensorflow-gpu \
    pytorch \
    torchvision \
    -c pytorch \
    && mamba clean -tipsy \
    && fix-permissions $CONDA_DIR \
    && fix-permissions /home/$NB_USER