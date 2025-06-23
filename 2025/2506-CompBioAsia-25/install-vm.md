# Software Installation for Compbio Asia 2025

Table of Content:
- [Required software](#required-sofware)
- [NOTES FOR OPERATORS](#notes-for-operators)
- [Installation](#Installation)
    - [Base](#base)
    - [CPU-only](#cpu-only)
    - [GPU](#gpu)
        - [AMBER24](#amber24)
        - [AmberTools25](#ambertools25)
        - [AutoDock Vina](#autodock-vina)
    - [Docker Images](#docker-images)
- [Data](#data)

---

## Required sofware

| Software | Version | Instructions/Documentation/Github | Notes |
|---|---|---|---|
| HMMER | 3.4 |	http://hmmer.org/documentation.html ||
| MMseqs2 | 17 | https://github.com/soedinglab/MMseqs2 ||
| nail | 0.3.0 | https://github.com/TravisWheelerLab/nail ||
| Pytorch | 2.7 | https://pytorch.org/ ||
| Python |3.11 | ||
| Gemmi || https://gemmi.readthedocs.io/en/latest/install.html#python-module ||
| Pandas ||||
| Matplotlib ||||
| Numpy	||||
| AlphaFold2 |||
| AlphaFold3 | | https://github.com/google-deepmind/alphafold3 ||
| AMBERTools25 | 25 | 	https://ambermd.org/GetAmber.php#amber ||
| AMBER24 | 24 | https://ambermd.org/GetAmber.php#amber ||
| AutoDock Vina ||| OPENCL installation due to infrastructure constraints |

> [!NOTE]
> Most of tools & dependencies (e.g., **conda**, **CUDA**) are installed in **`/opt/tools`**. **apt install** software is installed in default location. **Amber**, **Ambertools**, **AutoDock Vina** are installed in **`/data/tools`**.
> This is because:
> - **conda**: users might want to install further tools that may impact other attendees.
> - **CUDA**: requires to be "closer" to the VM and may cause issues if share isn't available.

---

## NOTES FOR OPERATORS

These notes are going to be in reverse chronological order (i.e., latest date comes first, oldest date comes last). These notes are expected to grow overtime. For installation notes, skip to the [Installation](#installation) section.

- **2025-06-20**:
    - 

- **2025-06-18**:
    - JupyterHub notes on attaching Manila share POST-deployment: [mounting_existing_manila_share.md at GitLab](https://gitlab.com/stack0/cacao-tf-jupyterhub/-/blob/main/docs/mounting_existing_manila_share.md?ref_type=heads#mounting-the-manila-share-into-jupyter-notebook). Kudos to Edwin!

- **2025-06-16**:
    - UPDATED SNAPSHOTS ARE: 
        - GPU: **`compbio-gpu-01-250616`**
        - CPU: **`compbio-cpu-01-250616`**
    - Made changes to CPU and GPU machines to address group ownership:
        - A new group is created: `sudo groupadd compbio-shared` (ran once in GPU and CPU new snapshots **DOES NOT REQUIRE TO BE EXECUTED AGAIN**)
        - Set ownership and permissions on `/data/user_dirs` (**DOES NOT REQUIRE TO BE EXECUTED AGAIN**):
            ```
            sudo chown -R root:compbio-shared /data/user_dirs       # Sets group ownership recursively
            sudo chmod -R g+rwX /data/user_dirs                     # Adds group read/write for files and execute on directories
            sudo find /data/user_dirs -type d -exec chmod g+s {} +  # Enables setgid bit on all directories so new files inherit group
            ```
        - Created `/etc/profile.d/add_to_compbio_shared.sh` that automatically runs at login (and made executable with `sudo chmod +x /etc/profile.d/add_to_compbio_shared.sh`):
            ```
            #!/bin/bash

            # Only for normal users, skip root and system users (UID<1000)
            if [ "$UID" -ge 1000 ] && [ "$EUID" -ne 0 ]; then

                # If not already in compbio-shared group
                if ! id -nG "$USER" | grep -qw compbio-shared; then
                    sudo /usr/local/bin/add-to-compbio-shared "$USER"
                fi

                # Set umask for group rw permissions on new files/dirs
                umask 002
            fi
            ```
        - Created a helper script `/usr/local/bin/add-to-compbio-shared` (and made executable with `sudo chmod +x /usr/local/bin/add-to-compbio-shared`):
            ```
            #!/bin/bash
            USER_TO_ADD="$1"
            if [ -z "$USER_TO_ADD" ]; then
                echo "Usage: $0 username"
                exit 1
            fi
            if id -nG "$USER_TO_ADD" | grep -qw compbio-shared; then
                exit 0
            fi
            usermod -aG compbio-shared "$USER_TO_ADD"
            ```
            - Allow passwordless sudo for the helper script: created file `/etc/sudoers.d/add-to-compbio-shared` (added right permissions with `sudo chmod 440 /etc/sudoers.d/add-to-compbio-shared`):
                ```
                ALL ALL=NOPASSWD: /usr/local/bin/add-to-compbio-shared
                ```
    - Made small script to delete clutter upon login (`/etc/profile.d/cleanup_home_dirs.sh`, and made executable `sudo chmod +x /etc/profile.d/cleanup_home_dirs.sh`):
        ```
        #!/bin/bash

        # Only apply to real user sessions (UID >= 1000, non-root)
        if [ "$UID" -ge 1000 ] && [ "$EUID" -ne 0 ]; then

        TARGET_DIRS=(Desktop Documents Downloads Music Pictures Public Templates Videos)

        found=false

        for dir in "${TARGET_DIRS[@]}"; do
            fullpath="$HOME/$dir"
            if [ -d "$fullpath" ]; then
            rm -rf "$fullpath"
            echo "Deleted: $fullpath"
            found=true
            fi
        done

        if ! $found; then
            echo "Nothing to delete in home"
        fi

        fi
        ```

- **2025-06-15**:
    - Added Docker section to this file.

- **2025-06-14**:
    - All VMs requested have been created.
    - **Base** image: **`compbio-base-01-250613`** (Ubuntu 22, 37G available disk space)
        - Contains the following software:
            - Conda
                - Python=3.11
                - Matplotlib
                - Numpy
                - Gemmi
                - Pytorch (GPU build)
            - tmux
            - nvtop
    - **CPU-only** image: **`compbio-cpu-00-250612`** (Ubuntu 24, 34G available disk space) 
        - Contains:
            - Everything from the **Base** installation (note Ubuntu 22 > 24, should not affect operations)
            - rustup
            - nail
            - HMMER
            - mmseqs2
    - **GPU** image: **`compbio-gpu-00-250614`** (Ubuntu 22, 27G available disk space)
        - Contains:
            - Everything from the **Base** installation
            - Conda
                - scipy
            - AMBER24
            - AmberTools25
            - AutoDock Vina (OpenCL) (executed with `autodock_gpu_128wi`)
    - Docker Images:
        - **Protein Ligand Prep**:
            - refer to the **`install-notes-jh-docker.md`** document.
            - Pull on a GPU capable machine: `cosimichele/jupyter-amber:gpu-250609`
        - Alphafold2: TBD
        - Alphafold3: TBD
    - **Data folder** (4TB):
        - Data can be deposited in `/data/user_dirs`.

---

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

**AMBER24**, **Ambertools25**, **AutoDock Vina** have their own subsections as both tools require extra care.

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
5. Create file that allow for Amber to be executed for everyone: `sudo nano /etc/profile.d/amber.sh`
    - Add: 
    ```
    source /data/tools/pmemd24/amber.sh
    export AMBERHOME=/data/tools/pmemd24/
    ```
    - Source: `source /etc/profile.d/amber.sh`
    - Doing `pmemd.cuda` runs.
    ```

      Unit  115 Error on OPEN: mdin                                                                                                                                                                                                                                                         
    STOP PMEMD Terminated Abnormally!
    ```

#### AmberTools25

>[!IMPORTANT]
> **Critical**: Follow the installation instructions in the [manual](https://ambermd.org/doc12/Amber25.pdf) found at https://ambermd.org/GetAmber.php

1. Manually Download **AmberTools25** (`ambertools25.tar.bz2`) from https://ambermd.org/GetAmber.php
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
4. Install **scipy**: `conda install scipy`
5. Execute the cmake install: `sudo ./run_cmake`
6. Install: `make install`
7. Add to `amber.sh` file that allow for Amber to be executed for everyone: `sudo nano /etc/profile.d/amber.sh`
    - Add: 
    ```
    source /data/tools/ambertools25/amber.sh
    export AMBERSOURCE=/data/tools/ambertools25_src
    ```
    The entire file should look like the following:
    ```
    source /data/tools/pmemd24/amber.sh
    export AMBERHOME=/data/tools/pmemd24/
    source /data/tools/ambertools25/amber.sh
    export AMBERSOURCE=/data/tools/ambertools25_src
    ```
    - Log out and log back in (close and reopen the VM)
    - Doing `pdb4amber` and `tleap` gives output

#### AutoDock Vina

AutoDock Vina with GPU support should be installed by compiling from source.

>[!IMPORTANT]
> Unfortunately, we are running a GRID A100 ("not a full GPU, only a part of it"; For the nerds: `GRID A100X-8C, based on the NVIDIA A100 architecture = SM_80`). AutoDock Vina fails if asked to be built for a full GPU (i.e., `make DEVICE=GPU TARGETS="80"`).
>
> Therefore we need to build using **OPENCL**.

>[!NOTE]
> AutoDock Vina is executed with **`autodock_gpu_128wi`** from the SHELL.

1. Clone the repository: `cd /data/raw_tools/ && git clone https://github.com/ccsb-scripps/AutoDock-GPU.git && cp -r /data/raw_tools/AutoDock-GPU /data/tools && cd /data/tools`
2. Export CUDA include and libraries to PATH: `export GPU_INCLUDE_PATH=/usr/local/cuda/include && export GPU_LIBRARY_PATH=/usr/local/cuda/lib64`
3. Run make: `make DEVICE=OPENCL TARGETS="80"`
4. Create file that allow for AutoDock to be executed for everyone: `sudo nano /etc/profile.d/autodock.sh`
    - Add: 
    ```
    export  PATH=$PATH:/data/tools/AutoDock-GPU/bin
    ```
    - Log out and log back in (close and reopen the VM)
5. Execute for testing:
    ```
    autodock_gpu_128wi --ffile /data/tools/AutoDock-GPU/input/1stp/derived/1stp_protein.maps.fld --lfile /data/tools/AutoDock-GPU/input/1stp/derived/1stp_ligand.pdbqt
    
    AutoDock-GPU version: v1.6-7-ga46ab564d2ac5f1a1523f65239b43505a1c29364

    Running 1 docking calculation

    Kernel source used for development:      ./device/calcenergy.cl
    Kernel string used for building:         ./host/inc/stringify.h
    Kernel compilation flags:                 -I ./device -I ./common -DN128WI   -cl-mad-enable
    OpenCL device:                           GRID A100X-8C
    (Thread 1 is setting up Job #1)

    Running Job #1
        Using heuristics: (capped) number of evaluations set to 1132076
        Local-search chosen method is: ADADELTA (ad)

    Executing docking runs, stopping automatically after either reaching 0.15 kcal/mol standard deviation of
    the best molecules of the last 4 * 5 generations, 42000 generations, or 1132076 evaluations:

    Generations |  Evaluations |     Threshold    |  Average energy of best 10%  | Samples | Best Inter + Intra
    ------------+--------------+------------------+------------------------------+---------+-------------------
            0 |          150 |  643.21 kcal/mol |  165.75 +/-  127.32 kcal/mol |       4 |   11.29 kcal/mol
            5 |        29635 |  643.21 kcal/mol |   18.99 +/-   58.86 kcal/mol |     149 |   -8.81 kcal/mol
            10 |        56368 |   24.92 kcal/mol |   -8.72 +/-    0.25 kcal/mol |      13 |   -9.21 kcal/mol
            15 |        82493 |   -8.43 kcal/mol |   -9.04 +/-    0.20 kcal/mol |       9 |   -9.37 kcal/mol
            20 |       108465 |   -8.70 kcal/mol |   -9.39 +/-    0.18 kcal/mol |       8 |   -9.61 kcal/mol
            25 |       134048 |   -9.05 kcal/mol |   -9.34 +/-    0.17 kcal/mol |      14 |   -9.63 kcal/mol
            30 |       159781 |   -9.16 kcal/mol |   -9.63 +/-    0.07 kcal/mol |       9 |   -9.72 kcal/mol
            35 |       185166 |   -9.38 kcal/mol |   -9.66 +/-    0.11 kcal/mol |      12 |   -9.92 kcal/mol
            40 |       210415 |   -9.52 kcal/mol |   -9.74 +/-    0.14 kcal/mol |      12 |   -9.97 kcal/mol
            45 |       235930 |   -9.56 kcal/mol |   -9.77 +/-    0.14 kcal/mol |       9 |  -10.00 kcal/mol
            50 |       261203 |   -9.54 kcal/mol |   -9.83 +/-    0.13 kcal/mol |      10 |  -10.02 kcal/mol
            55 |       286668 |   -9.64 kcal/mol |   -9.88 +/-    0.12 kcal/mol |      10 |  -10.03 kcal/mol
            60 |       312079 |   -9.70 kcal/mol |   -9.86 +/-    0.11 kcal/mol |      13 |  -10.03 kcal/mol
            65 |       337489 |   -9.74 kcal/mol |   -9.91 +/-    0.06 kcal/mol |      12 |  -10.03 kcal/mol
    ------------+--------------+------------------+------------------------------+---------+-------------------

                                    Finished evaluation after reaching
                                    -9.87 +/-    0.11 kcal/mol combined.
                                45 samples, best inter + intra   -10.03 kcal/mol.


    Job #1 took 0.229 sec after waiting 0.165 sec for setup

    (Thread 1 is processing Job #1)
    Run time of entire job set (1 file): 0.418 sec
    Processing time: 0.024 sec

    All jobs (1) ran without errors.
    ```

AutoDock Vina executes as requested.


#### Installation Complete

GPU image complete, snapshot created: `compbio-gpu-00-250614`.


### Docker Images

- **Amber MD Prep**: please refer to [Installation for Amber MD Prep](https://github.com/CosiMichele/workshop-materials/blob/main/2025/2506-CompBioAsia-25/amber_docker/install-notes-jh-amber-docker.md)(`install-notes-jh-amber-docker.md`).

### Alphafold 3

Alphafold3 has a number of requirements:

Hardware:
    - Access to  NVIDIA A100 80 GB GPU or  NVIDIA H100 80 GB GPU
    - 1TB disk space
    - 64GB RAM

Software:
    - CUDA 12.6
    - cuDNN
    - JAX
    - nvidia-ctk

Additionally, conda is installed (see steps above).

#### Installing requirements

1. **CUDA 12.6**:
    - Download keyring and update: `wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-keyring_1.1-1_all.deb && sudo dpkg -i cuda-keyring_1.1-1_all.deb && sudo apt-get update`
    - Install: `sudo apt-get -y install cuda-toolkit-12-6`
    - Add paths to `~/.bashrc`:
    ```
    export PATH=/usr/local/cuda-12.6/bin:$PATH
    export LD_LIBRARY_PATH=/usr/local/cuda-12.6/lib64:$LD_LIBRARY_PATH
    ```
    - `source ~/.bashrc`
    - Set default CUDA version by creating link
    ```
    sudo ln -sf /usr/local/cuda-12.6 /usr/local/cuda
    ```
    > [!NOTE]
    > Likely when doing `nvidia-smi` CUDA Version 12.4 is reported. However, when doing `nvcc -V` release 12.6 is showing.
    > ```
    > exouser@mentally-stirring-oyster:~$ nvidia-smi 
    > Sat Jun 21 15:06:35 2025       
    > +-----------------------------------------------------------------------------------------+
    > | NVIDIA-SMI 550.144.03             Driver Version: 550.144.03     CUDA Version: 12.4     |
    > |-----------------------------------------+------------------------+----------------------+
    > | GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
    > | Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
    > |                                         |                        |               MIG M. |
    > |=========================================+========================+======================|
    > |   0  NVIDIA H100 80GB HBM3          On  |   00000000:04:00.0 Off |                    0 |
    > | N/A   36C    P0             72W /  700W |      14MiB /  81559MiB |      0%      Default |
    > |                                         |                        |             Disabled |
    > +-----------------------------------------+------------------------+----------------------+
    >                                                                                         
    > +-----------------------------------------------------------------------------------------+
    > | Processes:                                                                              |
    > |  GPU   GI   CI        PID   Type   Process name                              GPU Memory |
    > |        ID   ID                                                               Usage      |
    > |=========================================================================================|
    > |    0   N/A  N/A      1398      G   /usr/lib/xorg/Xorg                              4MiB |
    > +-----------------------------------------------------------------------------------------+
    > exouser@mentally-stirring-oyster:~$ nvcc -V
    > nvcc: NVIDIA (R) Cuda compiler driver
    > Copyright (c) 2005-2024 NVIDIA Corporation
    > Built on Tue_Oct_29_23:50:19_PDT_2024
    > Cuda compilation tools, release 12.6, V12.6.85
    > Build cuda_12.6.r12.6/compiler.35059454_0
    > ```
    >
    > This mismatch is normal and expected behavior.
    > `nvidia-smi` reports the maximum CUDA version supported by the currently loaded NVIDIA driver, not the version of the CUDA toolkit installed.
    >
    > The driver version (550.144.03) is from the 12.4 driver series, and although it's compatible with CUDA 12.6, it does not declare that directly.
2. **cuDNN**:
    - Assuming that the previous steps have been taken, you can install with `sudo apt-get -y install cudnn` (or `sudo apt-get -y install cudnn-cuda-12` for specific version 12).
3. **JAX**:
    - Install using pip: `pip install --upgrade "jax[cuda12]"`
4. **nvidia-ctk** (following [these instructions](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md#installing-nvidia-support-for-docker)):
    - Download key: 
    ```
    curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg \
  && curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
    sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
    ```
    - Update and istall: `sudo apt-get update && sudo apt-get install -y nvidia-container-toolkit`
    - Set runtime: `nvidia-ctk runtime configure --runtime=docker --config=$HOME/.config/docker/daemon.json`
    - Restart docker: `systemctl --user restart docker`
    - Configure: `sudo nvidia-ctk config --set nvidia-container-cli.no-cgroups --in-place`


#### Obtaining AlfaFold3 Code and Data

From here onwards we can stick to the [official documentation](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md#obtaining-alphafold-3-source-code).

1. Clone repository: `git clone https://github.com/google-deepmind/alphafold3.git && cd alphafold3`
2. Obtain genetic databases
    > [!IMPORTANT]
    > **Make sure that you have enough disk space to run the script.**

---

## Data

Data will be stored in the `/data/` shared storage (4TB disk space). This is a [Manila](https://wiki.openstack.org/wiki/Manila)-type shared storage, where you can find:

- `/data/raw_tools`: raw files for installing tools. **DO NOT TOUCH**.
- `/data/raw_data`: location for raw data. **Do NOT TOUCH/WRITE**.
- `/data/tools`: installation location for Amber24, AmberTools25, AutoDock Vina. **DO NOT TOUCH**.
- `/data/user_dirs`: directory for users to create folders and save files.

---


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

1. install nvcc 
2. install cuDNN
3. install conda
4. install JAX
5. Download data
6. Build docker
7. Execute docker