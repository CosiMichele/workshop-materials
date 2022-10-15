#!/bin/bash

wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
chmod 775 Miniconda3-py39_4.12.0-Linux-x86_64.sh
./Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -u
~/miniconda3/bin/conda init
source ~/.bashrc
conda install -y mamba -n base -c conda-forge
mamba
