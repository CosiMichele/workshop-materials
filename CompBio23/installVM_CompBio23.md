# Installation and Instructions for all required packages

In order to address some of the concerns around scratch disk space, separate virtual machines (VMs) are set up addressing Genomics, MD, and Jupyter/RStudio requirements.

Access to these machines is on demand and/or according to the necessities of the day. Data will be stored in a shared directory accessible by all machines.

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
    -  Ambertools23 (https://ambermd.org/GetAmber.php) 
    -  Amber22/pmemd.cuda (this can wait – we’re figuring out the license)
    -  Gromacs2022.1 (https://manual.gromacs.org/current/download.html)
    -  g_mmpbsa (https://rashmikumari.github.io/g_mmpbsa/   - someone has suggested that this may come as part of Gromacs2022.1???)
- Screening
    - R (https://www.r-project.org/)
    -  PDB2PQR (https://sourceforge.net/projects/pdb2pqr/)
    -  AutoDock Vina (GPU version) (https://github.com/ccsb-scripps/AutoDock-Vina) (in drugsniffer)
    -  OpenBabel (http://openbabel.org/wiki/Main_Page)

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
        4. **Important**: because Tophat2 only runs on Python2, a conda environemnt had to be created for support. This was created with `conda create -n tophat-env python=2.7`. 
    - Test command: `conda activate tophat-env`, then `tophat --version`
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
    - **Alternative installation**:
        - `conda install -c bbglab gitools`
        - Does not work as intended (persisting Java issue).
- Cuda:
    - Version: `1.5`
    - Installation command(s):
        1. wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
        2. sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
        3. wget https://developer.download.nvidia.com/compute/cuda/11.5.0/local_installers/cuda-repo-ubuntu2004-11-5-local_11.5.0-495.29.05-1_amd64.deb
        4. sudo dpkg -i cuda-repo-ubuntu2004-11-5-local_11.5.0-495.29.05-1_amd64.deb
        5. sudo apt-key add /var/cuda-repo-ubuntu2004-11-5-local/7fa2af80.pub
        6. sudo apt-get update
        7. sudo apt-get -y install cuda
        8. export PATH=/usr/local/cuda-11.5/bin${PATH:+:${PATH}}
        9. export LD_LIBRARY_PATH=/usr/local/cuda-11.5/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
        10. sudo reboot
- Amber22, Ambertools22:
    - Version: Ambertools22
    - Installation command(s):
        1. Manually downloaded from https://ambermd.org/GetAmber.php (used own name/institution to download) + used licenced Amber22 (private downloa link)
        2. `bunzip2 AmberTools22.tar.bz2 && bunzip2 Amber22.tar.bz2`
        2. `tar -xvf AmberTools22.tar && tar -xvf Amber22.tar && rm Amber*`
        3. `cd amber22_src/build`
        4. Modify `./run_cmake` manually using editor such that `DCUDA=TRUE`, `-DDOWNLOAD_MINICONDA=FALSE` (miniconda already installed) **Important**: in case you want multithreading, please have `-DMPI=TRUE`
        6. `./run_cmake`
        7. `sudo make install`
        8. Added (with vim) `export AMBERSOURCE=/home/exouser/tools/amber22_master/amber22_src` and `source /home/exouser/tools/amber22_master/amber22/amber.sh` to `~/.bashrc`
    - Test command: `pdb4amber`, `antechamber`, `reduce`, `tleap` all give output.
- AmberTools23:
    - Version: `23`
    - Installation instructions: https://ambermd.org/doc12/Amber23.pdf
    - Installation procedure and commands:
        1. Manually downloaded source code from https://ambermd.org/GetAmber.php.
        2. Decompress with `tar xvfj AmberTools23.tar.bz2` (output is `amber22_src`). Change directory: `cd amber22_src/build`
        3. Modify `./run_cmake` manually using editor such that `DCUDA=TRUE`, `-DDOWNLOAD_MINICONDA=FALSE` (miniconda already installed) **Important**: in case you want multithreading, please have `-DMPI=TRUE`
        4. Build with `./run_cmake` and `sudo make install`
- GROMACS:
    - Version: `2022.1`
    - Installation command(s):
        1. `sudo apt-get  install -y libblas-dev liblapack-dev fftw3 fftw3-dev pkg-config`
        2. `wget ftp://ftp.gromacs.org/gromacs/gromacs-2022.1.tar.gz`
        3. `tar -xvf gromacs-2022.1.tar.gz && rm gromacs-2022.1.tar.gz`
        4. `mkdir build && cd build`
        5. `cmake .. -DGMX_EXTERNAL_BLAS=TRUE -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_FFT_LIBRARY=fftw3 -DGMX_MPI=on -DGMX_GPU=CUDA -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-11.5`
        6. `make`
        7. `make check`
        8. `sudo make install`
        9. `source /usr/local/gromacs/bin/GMXRC`
    - Test command: `gmx_mpi`
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
- DESeq2 (R):
    - Version: `DESeq2_1.36.0`
    - Installation command(s):
        1. `sudo apt-get -y install libcurl4-openssl-dev libxml2-dev`
        2. `sudo apt-get -y install libcurl4-gnutls-dev libxml2-dev libssl-dev`
        3. Open R as superuser `sudo R` (required to update packages)
        4. In R: `install.packages("BiocManager")`, yes to all options
        5. In R: `BiocManager::install("DESeq2")`
    - Test command: in R `library("DESeq2")` does not returns no Error.
- RStudio:
    - Version: `2022.02.3 Build 492`
    - Installation command(s): 
        1. Manual download from `https://www.rstudio.com/products/rstudio/download/#download`
        2. `sudo apt-get install -y libclang-dev`
        3. `sudo apt --fix-broken install -y`
        4. `sudo dpkg -i rstudio-2022.02.3-492-amd64.deb`
- Nextflow:
    - Version: `22.04.3.5703`
    - Installation command(s): 
        1. `curl -s https://get.nextflow.io | bash`
        2. `export PATH=$PATH:$(pwd)`
        3. Manually added `export PATH=$PATH:/home/exouser/tools` to `~/.bashrc`
    - Test command: `nextflow -v` 
    - **Important**: make sure you are in a folder you can write (e.g. Desktop)
- PDB2PQR:
    - Version: `3.5.2`
    - Installation command(s): 
        1. `git clone https://github.com/Electrostatics/pdb2pqr.git`
        2. `cd pbd2pqr && pip install .`
        3. `pip install pdb2pqr` 
    - **Important**: none of the releases are downloadable from the [link given](https://sourceforge.net/projects/pdb2pqr/). Tried installation following insructions from [readthedocs](https://pdb2pqr.readthedocs.io/en/latest/getting.html), unsure if installed correctly.
- AutoDock Vina:
    - Version: `1.2.3`
    - Installation command(s): Manually downloaded release from `https://github.com/ccsb-scripps/AutoDock-Vina`
    - `vina_split` installation: obtained `vina_split` from vina [1.1.2 release](https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz)
    - Test command: `vina_1.2.3_linux_x86_64 --version`
- OpenBabel:
    - Version: `3.1.0` (command line), `3.0.0` (GUI)
    - Installation command(s):
        1. Manually downloaded and decompressed latest release from `https://github.com/openbabel/openbabel/releases/tag/openbabel-3-1-1`
        2. `cd openbabel-3.1.1 && mkdir build && cd build`
        3. `cmake ..`
        4. `make -j2`
        5. `make test`
        6. `sudo make install`
        7. `sudo apt install openbabel-gui`
    - Test command: `obgui` opens the gui as expected, `obabel` runs in the command line

- RDKit:
    - Version: `2022.03.2`
    - Installation command(s):
        - **Important**: tried and failed to build from source, failed. Reverted to conda.
        1. `wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh`
        2. `chmod 775 Miniconda3-py39_4.12.0-Linux-x86_64.sh`
        3. `./Miniconda3-py39_4.12.0-Linux-x86_64.sh`
        4. `conda init`
        5. `conda create -c conda-forge -n my-rdkit-env rdkit`
    - Test commands:
        - `conda activate my-rdkit-env`
        - `python`
        - `from rdkit import Chem`
    - **Important**: installation for main occoured in instructor folder. Students require to repeat last step. 

---

## Software changelong between 2022 and 2023

- Removed: VMD, Chimera/ChimeraX, Autodock Vina, Drugsniffer, AlphaFold, FPADMET
- Updated: Amber22 -> Amber23, 
- Added: pldock (docker), MDTraj, PyMol

## Changelog (in slightly more detail)

- Removed docker images (Drugsniffer related)
- Removed from `/usr/local/tools`:
    - drug-sniffer
    - fpadmet
    - alphafold
- Renamed `/etc/profile.d/compbio22_init.sh` to `/etc/profile.d/compbio23_init.sh`
- Changes in `compbio23_init.sh`:
    - Removed:
        - `export FPADMET=/usr/local/tools/fpadmet/`
    - Added: 
        - `cp /etc/skel/.bashrc ~` (line 2)
        - `echo 'source /usr/local/tools/amber22_master/amber22/amber.sh' >> ~/.bashrc`
        - `conda init` (penultimate line)
        - `source ~/.bashrc` (last line)
- Enabled a smoother conda experience
    - Updated conda to a newer version (23.3.1)
    - installed mamba (`conda install -c conda-forge mamba`)