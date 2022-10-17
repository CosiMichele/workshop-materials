#!/bin/bash

cd ~/
#wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
#chmod +x Miniconda3-py39_4.12.0-Linux-x86_64.sh
#./Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -u
#~/miniconda3/bin/conda init
#source ~/.bashrc
conda install -y mamba -n base -c conda-forge
mamba

cd ~/
git clone https://github.com/annaprotasio/TE_ManAnnot
cd TE_ManAnnot && wget https://raw.githubusercontent.com/CosiMichele/workshop-materials/main/TE_symposium/te_annot_ubuntu.yml
conda env create -f te_annot_ubuntu.yml
chmod +x te_annot_ubuntu.yml
cd .. && mkdir Pfam_db && cd Pfam_db && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt
gunzip Pfam-A.hmm.gz && gunzip Pfam-A.hmm.dat.gz

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

head -n -4 ~/.bashrc > ~/.bashrc
