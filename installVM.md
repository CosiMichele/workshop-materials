# Installation and Instructions for all required packages

## Required software

- Genomics
    -  Bash and Zsh
    -  Python
    -  Java 8
    -  Autotools + make + cmake + gcc
    -  Blast (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
    -  HMMER (http://eddylab.org/software/hmmer/hmmer.tar.gz)
    -  MMSEQS2 (https://github.com/soedinglab/mmseqs2)
    -  MAFFT (https://mafft.cbrc.jp/alignment/software/)
    -  FastTree (http://www.microbesonline.org/fasttree/)
    -  PAML (http://abacus.gene.ucl.ac.uk/software/paml.html)
    -  Bowtie2 (https://github.com/BenLangmead/bowtie2)
    -  TopHat2 (https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz)
    -  Entrez Direct (https://www.ncbi.nlm.nih.gov/books/NBK179288/)
    -  SRA toolkit (https://hpc.nih.gov/apps/sratoolkit.html)
    -  Gitools (http://www.gitools.org/download)
- MD
    -  Bash and Zsh
    -  Ambertools22 (https://ambermd.org/AmberTools.php - no license required)
    -  Amber22/pmemd.cuda (this can wait – we’re figuring out the license)
    -  VMD (https://www.ks.uiuc.edu/Research/vmd/)
    -  ChimeraX (https://www.cgl.ucsf.edu/chimerax/)
    -  Gromacs2022.1 (https://manual.gromacs.org/current/download.html)
    -  g_mmpbsa (https://rashmikumari.github.io/g_mmpbsa/   - someone has suggested that this may come as part of Gromacs2022.1???)
- Screening
    - R (https://www.r-project.org/)
    - Drugsniffer (http://drugsniffer.org/)
    -  PDB2PQR (https://sourceforge.net/projects/pdb2pqr/)
    -  AutoDock Vina (GPU version) (https://github.com/ccsb-scripps/AutoDock-Vina) (in drugsniffer)
    -  OpenBabel (http://openbabel.org/wiki/Main_Page)
    -  ChimeraX (https://www.cgl.ucsf.edu/chimerax/)
    -  VMD (https://www.ks.uiuc.edu/Research/vmd/)
    -  AlphaFold (https://github.com/deepmind/alphafold)
    -  FPADMET - https://gitlab.com/vishsoft/fpadmet  (in drugsniffer)
    -  RDKIT (this requires Conda, I think) (https://www.rdkit.org/)   (in drugsniffer)

## Installation of software

Using a VM on JetStream2 to install software (`g3.small`), logged in via ssh.

- Base: Ubuntu 20.
- Prereq command execution:
    - `sudo apt-get update && sudo apt-get upgrade`
    - `mkdir ~/tools`

- Python:
    - Version: `v3.8.10`
    - Installation command(s): none, pre-installed with Base image.
    - Test command: `python 3 --version`
- Java 8:
    - Version: `openjdk 11.0.15 2022-04-19`
    - Installation command(s): `sudo apt-get install openjdk-8-jre`
    - Test command: `java --version`
- BLAST:
    - Version: `blast 2.13.0`
    - Installation command(s):
        1. `cd ~/tools`
        2. `wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz`
        3. `tar -xvf ncbi-blast-2.13.0+-x64-linux.tar.gz && rm ncbi-blast-2.13.0+-x64-linux.tar.gz`
        4. `export PATH=$PATH:$HOME/tools/ncbi-blast-2.13.0+/bin`
    - Test command: `blastn -version`
- HMMER:
    - Version: `HMMER 3.3.2`
    - Installation command(s):
        1. `wget http://eddylab.org/software/hmmer/hmmer.tar.gz`
        2. `tar -xvf hmmer.tar.gz && rm hmmer.tar.gz`
        3. `cd hmmr-3.3.2`
        4. `./configure --prefix /home/exouser/tools/hmmer-3.3.2`
        5. `make`
        6. `make check`
        7. `make install`
        8. `cd easel && make install`
        9. `export /home/exouser/tools/hmmer-3.3.2/bin`
    - Test command: `hmmscan -h`
    
- Tool:
    - Version:
    - Installation command(s):
    - Test command: