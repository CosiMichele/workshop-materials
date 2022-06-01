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

!!! Warning "Important"
        `PATH` is manually unpdated in `/etc/environment` whenever necessary


- Python2:
    - Version: `v2.7.18`
    - Installation command(s): `sudo apt install python`
    - Test command: `python --version`
- Python3:
    - Version: `v3.8.10`
    - Installation command(s): none, pre-installed with Base image.
    - Test command: `python3 --version`
- Java 8:
    - Version: `openjdk version "1.8.0_312"`
    - Installation command(s): 
        1. `sudo apt-get install openjdk-8-jdk`
        2. `sudo update-alternatives --config java` (selected 2: java 8)
    - Test command: `java -version`
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
- MMSEQS2:
    - Version: `96b2009982ce686e0b78e226c75c59fd286ba450`
    - Installation command(s):
        1. `wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz`
        2. `tar -xvf mmseqs-linux-avx2.tar.gz && rm mmseqs-linux-avx2.tar.gz`
        3. `export PATH=$(pwd)/mmseqs/bin/:$PATH`
    - Test command: `mmseqs --version`
- MAFFT:
    - Version: `v7.490`
    - Installation command(s):
        1. `wget https://mafft.cbrc.jp/alignment/software/mafft_7.490-1_amd64.deb`
        2. `dpkg -i mafft_7.490-1_amd64.deb && rm mafft_7.490-1_amd64.deb`
    - Test command: `mafft -v`
- FastTree:
    - Version: `2.1.11 Double precision (No SSE3)`
    - Installation command(s): `sudo apt-get install fasttree`
    - Test command: `FastTree`
- PAML:
    - Version: `4.9j`
    - Installation command(s): `sudo apt-get install paml`
    - Test command: `baseml`
- Bowtie 2:
    - Version: `2.3.5.1`
    - Installation command(s): `sudo apt install bowtie2`
    - Test command: `bowtie2 --version`
- Samtools:
    - Version: `1.10`
    - Installation command(s): `sudo apt-get install samtools`
    - Test command: `samtools --version`
- Tophat2:
    - Version: `v2.1.1`
    - Installation command(s):
        1. Manually downloaded `tophat-2.1.1.Linux_x86_64.tar.gz` using GUI from https://ccb.jhu.edu/software/tophat/downloads/ (`wget` was unable to establish SSL connection)
        2. `tar -xvf tophat-2.1.1.Linux_x86_64.tar.gz && rm tophat-2.1.1.Linux_x86_64.tar.gz`
        3. `export PATH=$(pwd)/tophat-2.1.1.Linux_x86_64:$PATH`
    - Test command: `tophat --version`
- Entrez Direct:
    - Installation command(s): 
        1. `sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"`
        2. `export PATH=${PATH}:${HOME}/edirect`
    - Test command: `esearch` (outputs `ERROR:  Missing -db argument`, which highlights esearch is installed)
- SRA-Toolkit:
    - Version: `3.0.0`
    - Installation command(s):
        1. `wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz`
        2. `tar -xvf sratoolkit.3.0.0-ubuntu64.tar.gz && rm sratoolkit.3.0.0-ubuntu64.tar.gz`
        3. `export PATH=$(pwd)/sratoolkit.3.0.0-ubuntu64/bin:$PATH`
        4. `vdb-config --interactive` (saved defaults and exited)
    - Test command: `fasterq-dump --version`
- GiTools:
    - Version: `2.3.1`
    - Installation command(s):
        1. `wget http://www.gitools.org/downloads/gitools-2.3.1-bin.zip`
        2. Unzipped with GUI
        3. `rm gitools-2.3.1-bin.zip` 
    - Test command: `./gitools`
    - **Important**: GiTools only accepts JAVA version 7. Although the [website](http://www.gitools.org/download) says that it works with Java 8, this isn't true. Furthermore, Java 7 is not supported on Ubuntu 16+. Further explainations [here](https://stackoverflow.com/questions/16263556/installing-java-7-on-ubuntu) and [here](https://askubuntu.com/questions/761127/how-do-i-install-openjdk-7-on-ubuntu-16-04-or-higher). 
- Ambertools22:
    - Version: Ambertools22
    - Installation command(s):
        1. Manually downloaded from https://ambermd.org/GetAmber.php (used own name/institution to download)
        2. `tar -xvf AmberTools22.tar.bz2 && rm AmberTools22.tar.bz2`
        3. `cd amber22_src/build`
        4. `./run_cmake`
        5. `make install`
        6. Added (with vim) `export AMBERSOURCE=/home/exouser/tools/amber22_src` and `source /home/exouser/tools/amber22/amber.sh` to `~/.bashrc`
    - Test command: `pdb4amber`, `antechamber`, `reduce`, `tleap` all give output. 
- VMD:
    - **Important**: more information required:
        - Download link: https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD
        - What version? 1.9.3 or 1.9.4 (alpha)
        - `Text mode` vs `text mode w/EGL` vs `OpenGL, CUDA, OptiX, OSPray`
- ChimeraX:
    - Version: `1.3`
    - Installation command(s):
        1. Manually downloaded from `https://www.cgl.ucsf.edu/chimerax/download.html` (accepted terms and conditions)
        2. `sudo apt-get install -y libxcb-xinerama0`
        3. `sudo apt --fix-broken install`
        4. `sudo dpkg -i ucsf-chimerax_1.3ubuntu20.04_amd64.deb && rm ucsf-chimerax_1.3ubuntu20.04_amd64.deb`
    - Test: Applciation opens when executed (from Desktop)
- GROMACS:
    - Version: `2022.1`
    - Installation command(s):
        1. `sudo apt-get install libblas-dev liblapack-dev fftw3 fftw3-dev pkg-config`
        2. `wget ftp://ftp.gromacs.org/gromacs/gromacs-2022.1.tar.gz`
        3. `tar -xvf gromacs-2022.1.tar.gz && rm gromacs-2022.1.tar.gz`
        4. `mkdir build && cd build`
        5. `cmake .. -DGMX_EXTERNAL_BLAS=TRUE -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_FFT_LIBRARY=fftw3`
        6. `make`
        7. `sudo make install`
        8. `source /usr/local/gromacs/bin/GMXRC`
    - Test command: `gmx -version`
- g_mmpbsa:
    - **Important**: more information required:
        - The latest version of GROMACS is 2022.1, whilst g_mmpbsa is only compatible with GROMACS <= 5.1. Thus a downgrade is necessary if this software is necessary.
- R:
    - Version: `4.2.0`
    - Installation command(s):
        1. `sudo apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common`
        2. `sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9`
        3. `sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'`
        4. `sudo apt install r-base`
    - Test command: `R --version`
- Tool:
    - Version:
    - Installation command(s):
    - Test command:
- Tool:
    - Version:
    - Installation command(s):
    - Test command:
- Tool:
    - Version:
    - Installation command(s):
    - Test command:
- Tool:
    - Version:
    - Installation command(s):
    - Test command:
- Tool:
    - Version:
    - Installation command(s):
    - Test command:
    
- Tool:
    - Version:
    - Installation command(s):
    - Test command: