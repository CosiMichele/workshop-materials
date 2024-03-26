# Install log for CompBio23

## Pre-setup
- Deployed a [g3.small](https://docs.jetstream-cloud.org/general/vmsizes/) VM through [CACAO](https://cacao.jetstream-cloud.org/), with Ubuntu22 as base
- Deployment has been given access to the compbio22 shared storage (`/data/`)
- Executed OS update & upgrade with `sudo apt-get update && sudo apt-get upgrade`

## Installation of software

- Conda/Mamba
    - Downloaded miniconda: `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`
    - Execute install: `./Miniconda3-latest-Linux-x86_64.sh`
    - Install location: `/data/tools/miniconda3`
    - Edited `/etc/environment`: added `/data/tools/miniconda3` to PATH
    - Sourced: `source ~/.bashrc`
    - Executed mamba install: `conda install -c conda-forge mamba`

- Java 8
    - Installation: `sudo apt-get install openjdk-8-jdk`
    - Selection of Java `sudo update-alternatives --config java` (selected 2: java 8)

- Cuda:
    - :exclamation::pencil: `nvidia-smi` and the `cuda` drivers are already installed but must be loaded first.
    - Load Cuda with: `module load nvhpc/22.11/nvhpc` (for versions of `NVIDIA-SMI 525.85.05`, `CUDA Version: 12.0` )
    - Check for `nvcc` with `nvcc --version` or `which nvcc`
    - Look at what other modules are available by doing `module avail`

- PyMol (Open source)
    - Intalled through conda/mamba: `mamba install -c conda-forge pymol-open-source`

- MDTraj
    - Installed with mamba: `mamba install -c conda-forge mdtraj`
    - Installed testing modules: `pip install pytest gsd`
    - Testing installation: `git clone https://github.com/mdtraj/mdtraj.git && cd mdtraj && py.test`

- BLAST
    - Download latest version of BLAST+: `wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.14.0+-x64-linux.tar.gz`
    - Decompress and remove compressed file: `tar -xvf ncbi-blast-2.14.0+-x64-linux.tar.gz && rm ncbi-blast-2.14.0+-x64-linux.tar.gz`
    - Export to path: `export PATH=$PATH:/data/tools/ncbi-blast-2.14.0+/bin`

- HMMER
    - Installable through conda/mamba: `mamba install -c bioconda hmmer`

- MMSEQ2
    - Installable through conda/mamba: `mamba install -c bioconda mmseq2`