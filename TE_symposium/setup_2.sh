#!/bin/bash

# Step 1: see https://github.com/CosiMichele/workshop-materials/blob/main/TE_symposium/setup_1.sh, lives in /etc/profile.d/

# Step 2: conda initiation and mamba installation
cd ~/
conda install -y mamba -n base -c conda-forge
mamba

# Step 3: Download TE_ManAnnot and related files
cd ~/
git clone https://github.com/annaprotasio/TE_ManAnnot
cd TE_ManAnnot && wget https://raw.githubusercontent.com/CosiMichele/workshop-materials/main/TE_symposium/te_annot_ubuntu.yml
conda env create -f te_annot_ubuntu.yml
chmod +x te_annot_ubuntu.yml
cd .. && mkdir Pfam_db && cd Pfam_db && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt
gunzip Pfam-A.hmm.gz && gunzip Pfam-A.hmm.dat.gz

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

# Step 6: hmm related commands
cd ~/Pfam_db
conda activate te_annot
hmmpress Pfam-A.hmm
conda deactivate
cd ~/
# one line alternative 
# cd ~/Pfam_db && conda activate te_annot && hmmpress Pfam-A.hmm & conda deactivate && cd ~/
