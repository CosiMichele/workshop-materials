Installed Conda
1. `wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh`
2. `chmod 775 Miniconda3-py39_4.12.0-Linux-x86_64.sh`
3. `./Miniconda3-py39_4.12.0-Linux-x86_64.sh`
4. `conda init`
5. `source ~/.bashrc`

DNAPipeTE Intallation:
1. sudo docker pull clemgoub/dnapipete:latest

TE_ManAnnot Installation:
1. `git clone https://github.com/annaprotasio/TE_ManAnnot`
2. cd TE_ManAnnot && conda create -f te_annot.yml 
