# Software Installation for Compbio Asia 2025

## Required sofware

| Software | Version | Instructions/Documentation/Github |
|---|---|---|
| HMMER | 3.4 |	http://hmmer.org/documentation.html |
| MMseqs2 | 17 | https://github.com/soedinglab/MMseqs2 |
| nail | 0.3.0 | https://github.com/TravisWheelerLab/nail |
| Pytorch | 2.7 | https://pytorch.org/ |
| Python |3.11 | |
| Gemmi || https://gemmi.readthedocs.io/en/latest/install.html#python-module |
| Pandas |||
| Matplotlib |||
| Numpy	|||
| AlphaFold3 | | https://github.com/google-deepmind/alphafold3 |
| AMBERTools25 | 25 | 	https://ambermd.org/GetAmber.php#amber |
| AMBER24 | 24 | https://ambermd.org/GetAmber.php#amber |
| AutoDock Vina |||
| AlphaFold2 |||

> [!NOTE]
> Most of tools & dependencies (e.g., **conda**, **CUDA**) are installed in **`/opt/tools`**. **Amber** and **Ambertools** are installed in **`/data/tools`**.
> This is because:
> - **conda**: users might want to install further tools that may impact other attendees.
> - **CUDA**: requires to be "closer" to the VM and may cause issues if share isn't available.

## Installation

Installation will be broken down into 2 machines: CPU-only and GPU-requiring machines.

All VMs will have the following (see [Base](#base)):
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- Python (3.11)
- Pandas
- Matplotlib
- Numpy
- Gemmi
- Pytorch (note: requires GPU VM for building)

For CPU only VMs, these sofware will be installed on top of universal software (see [CPU-only](#cpu-only)):
- HMMER
- MMSeq2
- nail

For GPU VMs, these sofware will be installed on top of universal software (see [GPU](#gpu)):
- AMBERTools25
- AMBER24
- AutoDock Vina

TBD:
- Alphafold2 
- Alphafold3

### Base

Building with flavour `g3.small`, Ubuntu 24.

1. Update system: `sudo apt-get update && sudo apt-get upgrade`
2. Create tools directory: `sudo mkdir /opt/tools`
3. Make it so that everyone can write to `tools`: `sudo chmod -R 777 /opt/tools`
4. **Conda** installation:
    - Download miniconda: `sudo wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/tools/miniconda.sh`
    - Execute installation: `cd /opt/tools && sudo chmod 777 ./miniconda.sh && sudo ./miniconda.sh` 
        - Select installation location: `/opt/tools/miniconda3`
        - Selected "yes" for conda initialization upon startup
        - Set up conda startup for any user; Create `sudo nano /etc/profile.d/conda.sh` and add:
        ```
        #!/bin/bash
        # System-wide Conda initialization + auto-activate base

        export PATH="/opt/tools/miniconda3/bin:$PATH"

        # Initialize Conda
        __conda_setup="$('/opt/tools/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
        if [ $? -eq 0 ]; then
            eval "$__conda_setup"
        else
            if [ -f "/opt/tools/miniconda3/etc/profile.d/conda.sh" ]; then
                . "/opt/tools/miniconda3/etc/profile.d/conda.sh"
            fi
        fi
        unset __conda_setup

        # Explicitly activate the base environment
        conda init
        ```
        - Make executable `sudo chmod +x /etc/profile.d/conda.sh`
    - Clean up: `sudo rm /opt/tools/miniconda.sh`
    - Allow for anyone to write/install software (not Best Practice, but for this purpose _it's ok_): `sudo chmod -R 777 /opt/tools/miniconda3`
5. Install **Python (3.11)**: `conda install python=3.11`
6. Install [**Base**](#base) packages (non GPU): `conda install -c conda-forge pandas matplotlib numpy gemmi -y`
7. Install **Pytorch**: `pip install torch --index-url https://download.pytorch.org/whl/cu118`
8. Install **tmux**: `sudo apt install tmux`
9. Install **nvtop**: `sudo apt install nvtop`

Base image complete, snapshot created: `compbio-base-00-250612`.

### CPU-only

Building using snapshot `compbio-base-00-250612` and an `m3.medium` flavour.

1. Install **HMMER**: `sudo apt install hmmer`
2. Install **MMseq2**: `conda install -c conda-forge -c bioconda mmseqs2`
3. Install **nail**:
    - Install `rustup`: 
    ```
    sudo env CARGO_HOME=/opt/tools/rust/cargo RUSTUP_HOME=/opt/tools/rust/rustup \
     curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sudo env \
     CARGO_HOME=/opt/tools/rust/cargo RUSTUP_HOME=/opt/tools/rust/rustup sh -s -- -y
    ```
    - Add to system-wide PATH: `sudo nano /etc/profile.d/rust.sh` and add
    ```
    export CARGO_HOME=/opt/tools/rust/cargo
    export RUSTUP_HOME=/opt/tools/rust/rustup
    export PATH=$CARGO_HOME/bin:$PATH
    ```
    - Make executable: `sudo chmod +x /etc/profile.d/rust.sh` and source `source /etc/profile.d/rust.sh`
    - Make folder writable by anyone (again, _it's ok_): `sudo chmod 777 -R /opt/tools/rust`
    - Install **nail**: `cargo install nail`

CPU image complete, snapshot created: `compbio-cpu-00-250612`.

### GPU

Building using snapshot `compbio-base-01-250613` and an `g3.small` flavour.

Note: downloading Amber files may take some time depending on network availability.

**AMBER24** and **Ambertools25** have their own subsections as both tools require extra care.

These tools require both CUDA drivers and NVIDIA Toolkit.

1. Install **CUDA 12.2 + NVIDIA 535**
    - Download the CUDA 12.2 local installer for Ubuntu 22.04: `wget https://developer.download.nvidia.com/compute/cuda/12.2.0/local_installers/cuda-repo-ubuntu2204-12-2-local_12.2.0-535.54.03-1_amd64.deb`
    - Install the repo and key:
        - Install deb package: `sudo dpkg -i cuda-repo-ubuntu2204-12-2-local_12.2.0-535.54.03-1_amd64.deb`
        - Install keyring: `sudo cp /var/cuda-repo-ubuntu2204-12-2-local/cuda-*-keyring.gpg /usr/share/keyrings/`
        - Update system: `sudo apt-get update`
    - Install CUDA: `sudo apt-get -y install cuda`
2. Create shell file so that `export` command is global:
    - Create file: `sudo nano /etc/profile.d/cuda.sh` and add:
    ```
    export PATH=/usr/local/cuda/bin:$PATH
    export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
    ```
    - Make executable: `sudo chmod +x /etc/profile.d/cuda.sh`
3. Cleanup deb: `rm cuda-repo-ubuntu2204-12-2-local_12.2.0-535.54.03-1_amd64.deb`

Logout and log back in, test with `nvidia-smi` and `nvcc --version`.

#### AMBER24

>[!IMPORTANT]
> **Critical**: Follow the installation instructions in the [manual](https://ambermd.org/doc12/Amber25.pdf) found at https://ambermd.org/GetAmber.php

1. Manually Download **AMBER24** (`pmemd24.tar.bz2`) from https://ambermd.org/GetAmber.php
    - This was done on a personal computer and transferred to the VM via `sftp` (`sftp put pmemd24.tar.bz2`)
2. Move `pmemd24.tar.bz2` to `/data/tools` with `sudo mv pmemd24.tar.bz2 /data/tools/`, and extract: `sudo bunzip2 pmemd24.tar.bz2 && sudo tar -xvf pmemd24.tar`. Results with `pmemd24_src` folder; `cd pmemd24_src`.
3. Open the `build` directory with `cd build` and edit `run_cmake` to reflect the following changes (lines 41-48):
    - `DCUDA=FALSE` -> `DCUDA=TRUE`
    - `DDOWNLOAD_MINICONDA=TRUE` -> `DDOWNLOAD_MINICONDA=FALSE`
    - Lines 41-48 should look like the following:
    ```
    cmake $AMBER_PREFIX/pmemd24_src \
        -DCMAKE_INSTALL_PREFIX=$AMBER_PREFIX/pmemd24 \
        -DCOMPILER=GNU  \
        -DMPI=FALSE -DCUDA=TRUE -DINSTALL_TESTS=TRUE \
        -DDOWNLOAD_MINICONDA=FALSE -DBUILD_PYTHON=FALSE \
        -DBUILD_PERL=FALSE -DBUILD_GUI=FALSE \
        -DPMEMD_ONLY=TRUE -DCHECK_UPDATES=FALSE \
        2>&1 | tee  cmake.log
    ```
    - This outputs the following:
    ```
    --                           3rd Party Libraries
    -- ---building bundled: -----------------------------------------------------
    -- blas - for fundamental linear algebra calculations
    -- lapack - for fundamental linear algebra calculations
    -- netcdf - for creating trajectory data files
    -- netcdf-fortran - for creating trajectory data files from Fortran
    -- kmmd - Machine-learning molecular dynamics
    -- ---using installed: ------------------------------------------------------
    -- zlib - for various compression and decompression tasks
    -- libbz2 - for various compression and decompression tasks
    -- libm - for fundamental math routines if they are not contained in the C library
    -- ---disabled: ------------------------------------------------
    -- libtorch - for fundamental math routines if they are not contained in the C library

    --                                Features:
    -- MPI:                               OFF
    -- MVAPICH2-GDR for GPU-GPU comm.:    OFF
    -- OpenMP:                            OFF
    -- CUDA:                              ON
    -- NCCL:                              OFF
    -- Build Shared Libraries:            ON
    -- Build GUI Interfaces:              OFF
    -- Build Python Programs:             OFF
    --  -Python Interpreter:              /opt/tools/miniconda3/bin/python (version 3.11)
    -- Build Perl Programs:               OFF
    -- Build configuration:               RELEASE
    -- Target Processor:                  x86_64
    -- Build Documentation:               ON
    -- Sander Variants:                   normal LES API LES-API QUICK-CUDA
    -- Install location:                  /data/tools/pmemd24/
    -- Installation of Tests:             ON

    --                               Compilers:
    --         C: GNU 11.4.0 (/usr/bin/gcc)
    --       CXX: GNU 11.4.0 (/usr/bin/g++)
    --   Fortran: GNU 11.4.0 (/usr/bin/gfortran)
    ```
4. Install: `make install`

#### Ambertools25

>[!IMPORTANT]
> **Critical**: Follow the installation instructions in the [manual](https://ambermd.org/doc12/Amber25.pdf) found at https://ambermd.org/GetAmber.php

1. Manually Download **Ambertools25** (`ambertools25.tar.bz2`) from https://ambermd.org/GetAmber.php
    - This was done on a personal computer and transferred to the VM via `sftp` (`sftp put ambertools25.tar.bz2`)
2. Move `ambertools25.tar.bz2` to `/opt/tools` with `sudo mv ambertools25.tar.bz2 /opt/tools/`, and extract: `sudo bunzip2 ambertools25.tar.bz2 && sudo tar -xvf ambertools25.tar`. Results with `ambertools25_src` folder; `cd ambertools25_src`.
3. Open the `build` directory with `cd build` and edit `run_cmake` to reflect the following changes (lines 40-45):
    - `DCUDA=FALSE` -> `DCUDA=TRUE`
    - `DDOWNLOAD_MINICONDA=FALSE` -> `DDOWNLOAD_MINICONDA=TRUE`
    - Lines 40-45 should look like the following:
    ```
    cmake $AMBER_PREFIX/ambertools25_src \
        -DCMAKE_INSTALL_PREFIX=$AMBER_PREFIX/ambertools25 \
        -DCOMPILER=GNU  \
        -DMPI=FALSE -DCUDA=TRUE -DINSTALL_TESTS=TRUE \
        -DDOWNLOAD_MINICONDA=FALSE \
        2>&1 | tee  cmake.log
    ```
4. Install nvidia-cuda-toolkit: `conda install -c "nvidia/label/cuda-12.2.0" cuda-toolkit`
6. Execute the cmake install: `sudo ./run_cmake`

