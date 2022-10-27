#!/bin/bash

# This script is in case the pfam server doesn't connect as needed.

# Step 4: Download Eccsplorer and related files
cd ~/
git clone https://github.com/crimBubble/ECCsplorer && cd ECCsplorer
mamba env create -f environment.yml
cd ~/miniconda3/envs/eccsplorer/bin
mv ~/ECCsplorer/ ECCsplorer
ln -s ECCsplorer/ECCsplorer.py eccsplorer
chmod +x eccsplorer
git clone https://bitbucket.org/petrnovak/repex_tarean.git && cd repex_tarean
conda activate eccsplorer
make && cd .. && ln -s repex_tarean/seqclust seqclust
chmod +x seqclust
conda deactivate

# Step 5: Download TETools and clean up
cd ~/
rm Miniconda3-py39_4.12.0-Linux-x86_64.sh*
git clone https://github.com/Dfam-consortium/TETools.git
sudo mv /etc/profile.d/uppsala2022_demoFiles ~/TETools/
sudo rm /etc/profile.d/setup_1.sh
head -n -4 ~/.bashrc > ~/.bashrc_tmp && mv ~/.bashrc_tmp ~/.bashrc
