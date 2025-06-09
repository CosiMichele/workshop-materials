FROM quay.io/jupyter/pytorch-notebook:cuda12-latest

USER root

# Install a few dependencies for text and monitoring instances
RUN apt update && \
    apt install -y lsb-release apt-transport-https curl libfreetype6-dev pkg-config libx11-dev gcc less software-properties-common apt-utils glances htop nano vim nvtop

RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - && \
    apt-get install -y nodejs

# Install and configure jupyter lab.
COPY jupyter_notebook_config.json /opt/conda/etc/jupyter/jupyter_notebook_config.json

# Add sudo to jovyan user
RUN apt update && \
    apt install -y sudo && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

ARG LOCAL_USER=jovyan
ARG PRIV_CMDS='/bin/ch*,/bin/cat,/bin/gunzip,/bin/tar,/bin/mkdir,/bin/ps,/bin/mv,/bin/cp,/usr/bin/apt*,/usr/bin/pip*,/bin/yum,/bin/snap,/bin/curl,/bin/tee,/opt'

RUN usermod -aG sudo jovyan && \
    echo "$LOCAL_USER ALL=NOPASSWD: $PRIV_CMDS" >> /etc/sudoers
RUN addgroup jovyan
RUN usermod -aG jovyan jovyan

# Update and clean package manager
RUN apt-get update --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER jovyan
WORKDIR /home/jovyan

RUN mamba

# Install required software
RUN mamba install -c conda-forge numpy cudatoolkit pdbfixer nglview openmm mdtraj ipywidgets=7.7.2 gitpython nbconvert nbformat requests packaging rich jax parmed netCDF4 openbabel widgetsnbextension jupyter_contrib_nbextensions nodejs -y
RUN pip install crossflow nglview[lab]
RUN nglview enable

# Rebuild the Jupyter Lab with new tools
RUN jupyter lab build

EXPOSE 8888

#COPY entry.sh /bin

#ENTRYPOINT ["bash", "/bin/entry.sh"]
