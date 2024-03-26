# Required software for TE Symposium:

## List of Software 

| Software | Installation/Test | Repository/Link |
| :---: | :---: | :---: |
| Miniconda & Mamba | [link](#conda-and-mamba) | [Miniconda](https://docs.conda.io/en/latest/miniconda.html) & [Mamba](https://github.com/mamba-org/mamba) |
| DNAPipeTE |[link](#dnapipete) | [repository](https://github.com/clemgoub/dnaPipeTE) |
| TE_ManAnnot | [link](#te_manannot) | [repository](https://github.com/annaprotasio/TE_ManAnnot) |
| ECCsplorer |[link](#eccsplorer) |[repository](https://github.com/crimBubble/ECCsplorer) |
| TETools | [link](#tetools) | [repository](https://github.com/Dfam-consortium/TETools) |
|||||

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

#### TE_ManAnnot:
- Installation:
    1. `git clone https://github.com/annaprotasio/TE_ManAnnot`
    2. `cd TE_ManAnnot && https://raw.githubusercontent.com/CosiMichele/workshop-materials/main/TE_symposium/te_annot_ubuntu.yml`
    3. `conda env create -f te_annot_ubuntu.yml`
    4. `cd .. && mkdir Pfam_db && cd Pfam_db && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt`
    5. `gunzip Pfam-A.hmm.gz && gunzip Pfam-A.hmm.dat.gz`
- Test:
    1. `hmmpress Pfam-A.hmm`

#### ECCsplorer:
- Installation:
    1. `git clone https://github.com/crimBubble/ECCsplorer && cd ECCsplorer`
    2. `mamba env create -f environment.yml`
    3. `cd ~/miniconda3/envs/eccsplorer/bin`
    4. `mv ~/ECCsplorer/ ECCsplorer`
    5. `ln -s ECCsplorer/ECCsplorer.py eccsplorer`
    6. `chmod +x eccsplorer`
    7. `git clone https://bitbucket.org/petrnovak/repex_tarean.git && cd repex_tarean`
    8. `conda activate eccsplorer`
    9. `make && cd .. &&  ln -s repex_tarean/seqclust seqclust`
    10. `chmod +x seqclust`
    11. `conda deactivate`
    12. sudo dpkg --add-architecture i386 && sudo apt-get update && sudo apt-get install libc6:i386 libncurses5:i386 libstdc++6:i386
- Test:
    1. `conda activate eccsplorer`
    2. `cd ~/miniconda3/envs/eccsplorer/bin/ECCsplorer/testdata`
    3. `efetch -db nucleotide -id CM009438.1 -seq_start 8971216 -seq_stop 9147030 -format fasta > chr1.fa && efetch -db nucleotide -id CM009440.1 -seq_start 44915437 -seq_stop 45118630 -format fasta > chr2.fa &&  efetch -db nucleotide -id CM009444.1 -seq_start 23920431 -seq_stop 24120625 -format fasta > chr3.fa`
    4. `cat chr1.fa chr2.fa chr3.fa | awk '/^>/{print ">chr" ++i; next}{print}' > RefGenomeSeq.fa`
    5. `rm chr1.fa chr2.fa chr3.fa`
    6. `efetch -db nucleotide -id JX455085.1 -format fasta > RefSeq_DB.fa`
    7. `cd ~/miniconda3/envs/eccsplorer/bin/ECCsplorer`
    8. `eccsplorer testdata/aDNA_R1.fastq testdata/aDNA_R2.fastq testdata/gDNA_R1.fastq testdata/gDNA_R2.fastq --reference_genome testdata/RefGenomeSeq.fa --database testdata/RefSeq_DB.fa --output_dir testrun --trim_reads tru3 --read_count 1000`

#### TETools:
- Installation:
    1. git clone https://github.com/Dfam-consortium/TETools.git
    2. docker pull dfam/tetools:1.5
    3. **Note**: dfam-tetools.sh looks for dfam/tetools:1.5. The VM is made to reflect this.

#### REPET:
- Installation:
    1. docker pull urgi/docker_vre_aio
