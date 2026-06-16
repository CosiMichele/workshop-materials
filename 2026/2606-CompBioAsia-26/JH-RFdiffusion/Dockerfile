FROM quay.io/jupyter/pytorch-notebook:cuda12-python-3.12.11

USER root

# Install a few dependencies for text and monitoring instances
RUN apt update && \
    apt install -y lsb-release apt-transport-https curl libfreetype6-dev pkg-config libx11-dev gcc less software-properties-common apt-utils glances htop nano vim nvtop unzip wget

RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
    apt-get install -y nodejs && \
    npm install -g npm@latest && \
    node --version && \
    npm --version

# Install and configure jupyter lab.
COPY jupyter_notebook_config.json /opt/conda/etc/jupyter/jupyter_notebook_config.json

# Install aria2
RUN apt-get install -y aria2

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

# Install Python packages
RUN pip install huggingface jedi omegaconf hydra-core icecream pyrsistent pynvml decorator numpy matplotlib ipywidgets py3Dmol transformers

# Install packages from git
RUN pip install git+https://github.com/NVIDIA/dllogger#egg=dllogger
RUN pip install --no-dependencies dgl -f https://data.dgl.ai/wheels/torch-2.4/cu124/repo.html
RUN pip install --no-dependencies e3nn==0.5.5 opt_einsum_fx
RUN pip install git+https://github.com/sokrypton/ColabDesign.git
RUN pip install git+https://github.com/Biohub/esm.git@main

# Clone RFdiffusion repository
RUN git clone https://github.com/sokrypton/RFdiffusion.git

# Install SE3Transformer
RUN cd RFdiffusion/env/SE3Transformer && pip install . && cd ../

# Download ananas tool
RUN wget -qnc https://files.ipd.uw.edu/krypton/ananas -O ananas && \
    chmod +x ananas

# Create params directory
RUN mkdir -p params

# Download RFdiffusion parameters
RUN mkdir -p RFdiffusion/models && \
    aria2c -q -x 16 http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt -o RFdiffusion/models/Base_ckpt.pt && \
    aria2c -q -x 16 http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt -o RFdiffusion/models/Complex_base_ckpt.pt && \
    aria2c -q -x 16 http://files.ipd.uw.edu/pub/RFdiffusion/f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt -o RFdiffusion/models/Complex_beta_ckpt.pt

# Download AlphaFold parameters
RUN cd params && \
    aria2c -q -x 16 \
    --allow-overwrite=true \
    --auto-file-renaming=false \
    https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar && \
    tar -xf alphafold_params_2022-12-06.tar && \
    rm alphafold_params_2022-12-06.tar && cd ~

# Download schedules
RUN aria2c -q -x 16 https://files.ipd.uw.edu/krypton/schedules.zip -o schedules.zip && \
    unzip -q schedules.zip && \
    rm schedules.zip

# Set environment variables
ENV DGLBACKEND=pytorch
ENV PYTHONPATH=/home/jovyan/RFdiffusion:$PYTHONPATH

# Rebuild the Jupyter Lab with new tools
RUN jupyter lab build

EXPOSE 8888

