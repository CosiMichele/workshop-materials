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
    - PDB2PQR (https://sourceforge.net/projects/pdb2pqr/)
    - AutoDock Vina (GPU version) (https://github.com/ccsb-scripps/AutoDock-Vina)
    - OpenBabel (http://openbabel.org/wiki/Main_Page)
    - RDKIT (this requires Conda, I think) (https://www.rdkit.org/)

## 2023 Additions

- AmberTools23:
    - Version: `23`
    - Installation instructions: https://ambermd.org/doc12/Amber23.pdf
    - Installation procedure and commands:
        1. Manually downloaded source code from https://ambermd.org/GetAmber.php.
        2. Decompress with `tar xvfj AmberTools23.tar.bz2` (output is `amber22_src`). Change directory: `cd amber22_src/build`
        3. Modify `./run_cmake` manually using editor such that `DCUDA=TRUE`, `-DDOWNLOAD_MINICONDA=FALSE` (miniconda already installed) :exclamation::pencil: **NOTE**: in case you want multithreading, please have `-DMPI=TRUE`.
        4. Install packages through conda: `mamba install -c conda-forge numpy scipy matplotlib`
        5. Build with `./run_cmake` and `sudo make install`
        6. Manually added `export AMBERSOURCE=/home/exouser/tools/amber22_master/amber22_src` and `echo 'source /usr/local/tools/amber22_master/amber22/amber.sh' >> ~/.bashrc` to `/etc/profile.d/compbio23_init.sh`. This will ensure that everytime the machine is turned on, Amber22 and AmberTools23 are executable from the path. 
        7. Made all directories and files in `$AMBERSOURCE` executable by all users with `sudo find ./* -type d -exec chmod 0777 {} \+` and  `sudo find ./* -type f -exec chmod 0777 {} \+` (this was necessary as otherwise the test would fail.)
    - Test:
        - `cd $AMBERSOURCE` and executed tests with `make test.serial` and `exportCUDA_VISIBLE_DEVICES=0 && make test.cuda.serial`. `make test.serial` executes with no errors, `test.cuda.serial` returns 2 erros, probably tied to GPU memory leaks (GPU out of memory)
- AutoDock Vina & Vina split:
    - Vina version: `vina_1.2.5_linux_x86_64`
    - Split version: `vina_split_1.2.5_linux_x86_64`
    - Installation procedure and commands: 
        1. Manually downloaded releases from `https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.5`
        2. Added path to `PATH` and renamed `vina_1.2.5_linux_x86_64` to `vina` and `vina_split_1.2.5_linux_x86_64` to `vina-split`
    - Test commands: `vina --version` and `vina-split --version` return version numbers and are executable from anywhere
    - :exclamation::pencil: **NOTE**: In case there's need to install the GPU version, go here: https://github.com/ccsb-scripps/AutoDock-GPU
- MDTraj:
    - Version: `1.9.8`
    - Installation instructions: https://www.mdtraj.org/1.9.8.dev0/installation.html#testing-your-installation
    - Intallation procedure and commands:
        1. `mamba install -c conda-forge mdtraj` (if mamba not installed use conda: `conda install -c conda-forge mdtraj`)
        2. `pip install pytest gsd`
    - Test:
        1. Executed `git clone https://github.com/mdtraj/mdtraj.git && cd mdtraj && py.test`. Test output resulted in 1000+ passed tests and 10 errors related to `rb` (read binary mode). Might need to look into this further.
- PyMOL:
    - Version: `2.5.0`
    - Source/repository: https://github.com/schrodinger/pymol-open-source
    - Installation instructions: https://pymolwiki.org/index.php/Linux_Install
    - Installation procedure and commands:
        1. Install requirements with `sudo apt-get install git build-essential python3-dev libglew-dev libpng-dev libfreetype6-dev libxml2-dev libmsgpack-dev python3-pyqt5.qtopengl libglm-dev libnetcdf-dev`
        2. Clone repositories: `git clone https://github.com/schrodinger/pymol-open-source.git`, `git clone https://github.com/rcsb/mmtf-cpp.git`, `mv mmtf-cpp/include/mmtf* pymol-open-source/include/`, `cd pymol-open-source`
        3. Set prefix location: `export prefix=/usr/local/tools/pymol-open-source/pymol-open-source-build`
        4. Install using the following: `python3 setup.py build install --home=$prefix`
    - Test:
        - `pymol-open-source-build/bin/pymol` results with an Error: `ImportError: /lib/x86_64-linux-gnu/libp11-kit.so.0: undefined symbol: ffi_type_pointer, version LIBFFI_BASE_7.0`. Looking more into it.
     
## Installation of software

Using a VM on JetStream2 to install software (`g3.small`), logged in via ssh.

- Base: Ubuntu 20.
- Prereq command execution:
    - `sudo apt-get update && sudo apt-get upgrade`
    - `mkdir /data/tools`

:warning: Warning "Important"
        `PATH` is manually unpdated in `/etc/environment` whenever necessary

- Python2:
    - Installed in own enviroment via conda
- Python3:
    - Version: `v3.10.10`
    - Installation command(s): none, installed with Conda.
    - Test command: `python3 --version`
- Java 8:
    - Version: `openjdk version "1.8.0_312"`
    - Installation command(s): 
        1. `sudo apt-get install openjdk-8-jdk`
        2. `sudo update-alternatives --config java` (selected 2: java 8)
    - Test command: `java -version`
- BLAST:
    - Version: `blast 2.14.0`
    - Installation command(s):
        1. `cd /data/tools`
        2. `wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.14.0+-x64-linux.tar.gz`
        3. `tar -xvf ncbi-blast-2.14.0+-x64-linux.tar.gz && rm ncbi-blast-2.14.0+-x64-linux.tar.gz`
        4. `export PATH=$PATH:/data/tools/ncbi-blast-2.14.0+/bin`
    - Test command: `blastn -version`
- HMMER:
    - Version: `HMMER 3.3.2`
    - Installation command(s):
        1. `mamba install -c bioconda hmmer`
    - Test command: `hmmscan -h`
- MMSEQS2:
    - Version: `13.45111`
    - Installation command(s):
        1. `mamba install -c bioconda mmseq2`
    - Test command: `mmseqs --version`
- MAFFT:
    - Version: `v7.520`
    - Installation command(s):
        1. `mamba install -c bioconda mafft`
    - Test command: `mafft -v`
- FastTree:
    - Version: `2.1.11 Double precision (No SSE3)`
    - Installation command(s): `mamba install -c bioconda fasttree`
    - Test command: `fasttree`
- PAML:
    - Version: `4.10`
    - Installation command(s): `mamba install -c bioconda paml`
    - Test command: `baseml`
- Bowtie 2:
    - Version: `2.5.1`
    - Installation command(s): `mamba install -c bioconda bowtie2`
    - Test command: `bowtie2 --version`
- Samtools:
    - Version: `1.13`
    - Installation command(s): `sudo apt-get install samtools`
    - Test command: `samtools --version`
- Tophat2:
    - Version: `v2.1.1`
    - Installation command(s):
        1. Manually downloaded `tophat-2.1.1.Linux_x86_64.tar.gz` using GUI from https://ccb.jhu.edu/software/tophat/downloads/ (`wget` was unable to establish SSL connection)
        2. `tar -xvf tophat-2.1.1.Linux_x86_64.tar.gz && rm tophat-2.1.1.Linux_x86_64.tar.gz`
        3. `exportex PATH=$(pwd)/tophat-2.1.1.Linux_x86_64:$PATH`
        4. **Important**: because Tophat2 only runs on Python2, a conda environemnt had to be created for support. This was created with `conda create -n tophat-env python=2.7`. 
    - Test command: `conda activate tophat-env`, then `tophat --version`
- Entrez Direct:
    - Version: `16.2`
    - Installation command(s): `mamba install -c bioconda entrez-direct`
    - Test command: `esearch -h`
- SRA-Toolkit:
    - Version: `3.0.5`
    - Installation command(s):`mamba install -c bioconda sra-tools`
    - Test command: `fasterq-dump --version`
- GiTools: ######### Followed alternative installation, still results in issues
    - Version: `2.3.1`
    - Installation command(s):
        1. `wget http://www.gitools.org/downloads/gitools-2.3.1-bin.zip`
        2. Unzipped with GUI
        3. `rm gitools-2.3.1-bin.zip`
    - Test command: `./gitools`
    - **Important**: GiTools only accepts JAVA version 7. Although the [website](http://www.gitools.org/download) says that it works with Java 8, this isn't true. Furthermore, Java 7 is not supported on Ubuntu 16+. Further explainations [here](https://stackoverflow.com/questions/16263556/installing-java-7-on-ubuntu) and [here](https://askubuntu.com/questions/761127/how-do-i-install-openjdk-7-on-ubuntu-16-04-or-higher). 
    - **Alternative installation**:
        - `conda create -n gitools openjdk=8`
        - `conda install -c bbglab gitools`
        - Does not work as intended (persisting Java issue).
- Cuda:
    - see [CompBio22 instructions](../CompBio22/installVM_CompBio22.md)
- Amber22:
    - Version: Amber22
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
    - :exclamation::pencil: **NOTE**: `pmemd.cuda` is installed and executable.
- GROMACS:
    - see [CompBio22 instructions](../CompBio22/installVM_CompBio22.md)
- Nextflow:
    - Version: `23.04.2.5870`
    - Installation command(s): 
        1. set Java to 11 with `sudo update-alternatives --config java` and select 0
        1. `wget -qO- https://get.nextflow.io | bash`
        2. `export PATH=$PATH:$(pwd)`
    - Test command: `nextflow -v` 
- PDB2PQR:
    - Version: `3.6.1`
    - Installation command(s): `mamba install -c conda-forge pdb2pqr`
    - Test command: pdb2pqr --version
- OpenBabel:
    - Version: `3.1.0` (command line), `3.0.0` (GUI)
    - Installation command(s):
        1. Manually downloaded and decompressed latest release from `https://github.com/openbabel/openbabel/releases/tag/openbabel-3-1-1`
        2. `cd openbabel-3.1.1 && mkdir build && cd build`
        3. `cmake .. -DCMAKE_INSTALL_PREFIX=/data/tools/openbabel -DBUILD_GUI=ON`
        4. `make -j2`
        5. `sudo make install`
        6. `export PATH=$PATH:/data/tools/openbabel/bin`
        7. `sudo apt install openbabel-gui`
    - Test command: `obgui` opens the gui as expected, `obabel` runs in the command line
- RDKit:
    - Version: `2022.03.2`
    - Installation command(s): `mamba install -c conda-forge rdkit`
    - Test commands:
        - `conda activate my-rdkit-env`
        - `python`
        - `from rdkit import Chem`

The following packages will be installed in the Dockerfile rather than the VM, although R and RStudio are available.
- R, RStudio:
    - R
        - Versions: R = `4.3.1` RStudio = `rstudio-2022.07.2+576 `
        - sudo apt update -qq
        - apt install --no-install-recommends software-properties-common dirmngr
        - wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
        - sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
        - sudo apt install --no-install-recommends r-base
    - RStudio
        - wget https://s3.amazonaws.com/rstudio-ide-build/electron/focal/amd64/rstudio-2023.06.1-524-amd64.deb
        - sudo dpkg -i rstudio-2023.06.1-524-amd64.deb
- DESeq2 (R):
    - Version: `DESeq2_1.36.0`
    - Installation command(s):
        1. `sudo apt-get -y install libcurl4-openssl-dev libxml2-dev`
        2. `sudo apt-get -y install libcurl4-gnutls-dev libxml2-dev libssl-dev`
        3. Open R as superuser `sudo R` (required to update packages)
        4. In R: `install.packages("BiocManager")`, yes to all options
        5. In R: `BiocManager::install(version = "3.17")`
        5. In R: `BiocManager::install("DESeq2")`
    - Test command: in R `library("DESeq2")` does not returns no Error.
- Bambu (R):
    - BiocManager::install("bambu")

---

## Software changelong

- Removed: 
    - VMD
    - Chimera/ChimeraX 
    - Drugsniffer 
    - AlphaFold 
    - FPADMET
    - g_mmpbsa
- Updated: 
    - AmberTools22 -> AmberTools23 
    - Autodock Vina 1.2.3 -> 1.2.5
- Added: 
    - pldock (docker) 
    - MDTraj
    - PyMol
- Issue:
    - Gromacs. Might want to use an older version of the VM to access.

### Changelog (in slightly more detail)

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

---

Note: [Edwin has good notes on shared storage for JH-vm4w](https://gitlab.com/stack0/jh-share-mounter)