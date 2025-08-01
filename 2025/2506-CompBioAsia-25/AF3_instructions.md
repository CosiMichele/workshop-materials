# AlphaFold3 Installation Instructions

Find the official [AF3 installation instructions](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md) in the [AF3 official repository](https://github.com/google-deepmind/alphafold3).

This guide takes into consideration that you have access to [JetStream2](https://docs.jetstream-cloud.org/)'s g5.xl flavour in order to install and execute AF3. Read more on JS2 flavours [here](https://docs.jetstream-cloud.org/general/instance-flavors/?h=flavo#jetstream2-gpu).

Table of Contents:
- [AlphaFold3 Installation Instructions](#alphafold3-installation-instructions)
  - [Installation Requirements](#installation-requirements)
  - [Installation](#installation)
    - [1. Install NVCC (NVIDIA CUDA Compiler)](#1-install-nvcc-nvidia-cuda-compiler)
    - [2. Install cuDNN \& JAX](#2-install-cudnn--jax)
    - [3. Clone Repository and Obtain Data](#3-clone-repository-and-obtain-data)
    - [4. Build Docker and Execute AF3](#4-build-docker-and-execute-af3)

## Installation Requirements

**Hardware**:
    - Access to  NVIDIA A100 80 GB GPU or NVIDIA H100 80 GB GPU
    - 1TB disk space
    - 64GB RAM

**Software**:
    - CUDA 12.6
      - Note: AF3 potentially *may* work with CUDA 12.4, but 12.6 is recommended.
    - cuDNN
    - JAX

**Notes on Models & Databases**:
    - [Model parameters](https://github.com/google-deepmind/alphafold3?tab=readme-ov-file#obtaining-model-parameters) can only be acquired by requested access. Please contact the AF3 team.
    - Databases require ~630GB of space; please ensure there is enough disk space to download and decompress data.

## Installation

These steps are done using an Ubuntu 22.04 distributiuon running on a g5.xl flavour JS2 VM. 

### 1. Install NVCC (NVIDIA CUDA Compiler)

1. Download keyring and update: `wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-keyring_1.1-1_all.deb && sudo dpkg -i cuda-keyring_1.1-1_all.deb && sudo apt-get update`
2. Install: `sudo apt-get -y install cuda-toolkit-12-6`
3. Add paths to `~/.bashrc`:
    ```
    export PATH=/usr/local/cuda-12.6/bin:$PATH
    export LD_LIBRARY_PATH=/usr/local/cuda-12.6/lib64:$LD_LIBRARY_PATH
    ```
4. `source ~/.bashrc`
5. Set default CUDA version by creating link
    ```
    sudo ln -sf /usr/local/cuda-12.6 /usr/local/cuda
    ```
6. Test by doing `nvcc --version`

>[!NOTE]
> There is a chance that when doing `nvidia-smi` the reported CUDA version is *NOT* 12.6. As long as `nvcc --version` reports version 12.6, AF3 will function a desired.


### 2. Install cuDNN & JAX

Assuming that the previous section ([Install NVCC](#1-install-nvcc-nvidia-cuda-compiler)) has taken place, there is no need to add the keyring. Otherwise, please carry out point 1 of the previous section.

Install cuDNN with `sudo apt-get -y install cudnn` 

>[!NOTE]
> One can specity the cuda version with`sudo apt-get -y install cudnn-cuda-12` for specific version 12 -- the standard installation should install the latest version.

Install [JAX](https://docs.jax.dev/en/latest/index.html) using pip: `pip install --upgrade "jax[cuda12]"`

### 3. Clone Repository and Obtain Data

This section takes into consideration a large directory (>1TB) named `data` and located in `/` : `/data`. From here onwards we can follow the [official documentation](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md#obtaining-alphafold-3-source-code).

1. Change directory to `/data`: `cd /data`
2. Clone repostory and change directory: `git clone https://github.com/google-deepmind/alphafold3.git && cd alphafold3`
3. Create folder for databases: `mkdir /data/public_databases`
4. Download databases: `./fetch_databases.sh /data/public_databases`
   - The command is `./fetch_databases.sh [<DB_DIR>]`, if `[<DB_DIR>]` is left empty, the databases are downloaded into `$HOME/public_databases`
   - Once the script has finished, you should have the following directory structure:
    ```
    mmcif_files/  # Directory containing ~200k PDB mmCIF files.
    bfd-first_non_consensus_sequences.fasta
    mgy_clusters_2022_05.fa
    nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta
    pdb_seqres_2022_09_28.fasta
    rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta
    rnacentral_active_seq_id_90_cov_80_linclust.fasta
    uniprot_all_2021_04.fa
    uniref90_2022_05.fa
    ```
5. Model paramenters:
    - As mentioned before, you will require to contact the AF3 team for the model parameters.
    - Once obtained, place model parameters in `/data/mod_params`.

### 4. Build Docker and Execute AF3

1. From within the `/data/alphafold3` directory, build Docker container with `docker build -t alphafold3 -f docker/Dockerfile .`
2. Create directories `/data/af_input` and `/data/af_output`; place [AF3 JSON](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md) inputs in `/data/af_input`.
3. Change directory to `/data/` (`cd /data/`) and execute with
    ```
    docker run -it \
        --volume ./af_input:/root/af_input \
        --volume ./af_output:/root/af_output \
        --volume ./mod_params:/root/models \
        --volume ./public_databases:/root/public_databases \
        --gpus all \
        alphafold3 \
        python run_alphafold.py \
        --json_path=/root/af_input/fold_input.json \
        --model_dir=/root/models \
        --output_dir=/root/af_output
    ```
