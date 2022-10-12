# Required software for TE Symposium:

## List of Software 

| Software | Installation/Test | Repository/Link |
| :---: | :---: | :---: |
| Miniconda & Mamba | [link](#conda-and-mamba) | [Miniconda](https://docs.conda.io/en/latest/miniconda.html) & [Mamba](https://github.com/mamba-org/mamba) |
| DNAPipeTE |[link](#dnapipete) | [repository](https://github.com/clemgoub/dnaPipeTE) |
| TE_ManAnnot | [link](#te_manannot) | [repository](https://github.com/annaprotasio/TE_ManAnnot) |

## Installation and testing

#### Conda and Mamba
1. `wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh`
2. `chmod 775 Miniconda3-py39_4.12.0-Linux-x86_64.sh`
3. `./Miniconda3-py39_4.12.0-Linux-x86_64.sh`
4. `conda init`
5. `source ~/.bashrc`
6. `conda install mamba -n base -c conda-forge`
7. `mamba`

#### DNAPipeTE:
- Installation:
    1. `docker pull clemgoub/dnapipete:latest`
- Test:
    1. `mkdir ~/Project && cd ~/Project`
    2. `docker run -it -v ~/Project:/mnt clemgoub/dnapipete:latest`
    3. `python3 dnaPipeTE.py -input /mnt/reads_input.fastq -output /mnt/output -RM_lib ../RepeatMasker/Libraries/RepeatMasker.lib -genome_size 170000000 -genome_coverage 0.1 -sample_number 2 -RM_t 0.2 -cpu 2`
    - **NOTE**: failed due to lack of test file.

#### TE_ManAnnot:
- Installation:
    1. `git clone https://github.com/annaprotasio/TE_ManAnnot`
    2. `cd TE_ManAnnot && wget https://github.com/CosiMichele/workshop-materials/blob/main/TE_symposium/te_annot_ubuntu.yml`
    3. `cd .. && mkdir Pfam_db && cd Pfam_db && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt`
    4. `gunzip Pfam-A.hmm.gz && gunzip Pfam-A.hmm.dat.gz`
- Test:
    1. `hmmpress Pfam-A.hmm`
