# Software Installation for Compbio 2025

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
2. Create tools directory: `mkdir /opt/tools`
3. Make it so that everyone can write to `tools`: `sudo chmod -R 777 /opt/tools`
4. **Conda** installation:
    - Download miniconda: ` wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/tools/miniconda.sh`
    - Execute installation: `cd /opt/tools && sudo chmod 777 ./miniconda.sh && sudo ./miniconda.sh` 
        - Select installation location: `/opt/tools/miniconda3`
        - Selected "yes" for conda initialization upon startup
        - Set up conda startup for any user; Create `/etc/profile.d/conda.sh` and add:
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
        - Add to `/etc/bash.bashrc` with `sudo vim /etc/bash.bashrc`:
        ```
        # >>> Force Conda base activation for all users >>>
        export PATH="/opt/tools/miniconda3/bin:$PATH"
        conda init
        # <<< End Conda base activation <<<
        ```
    - Clean up: `sudo rm /opt/tools/miniconda.sh`
    - Allow for anyone to write/install software (not Best Practice, but for this purpose _it's ok_): `sudo chmod -R 777 /opt/tools/miniconda3`
5. Install **Python (3.11)**: `conda install python=3.11`
6. Install [**Base**](#base) packages (non GPU): `conda install -c conda-forge pandas matplotlib numpy gemmi -y`
7. Install **`nvcc`**: `sudo apt install cuda-toolkit`
   - Note, this output:
   ```
   Restarting services...
    systemctl restart nvidia-persistenced.service
    Job for nvidia-persistenced.service failed because the control process exited with error code.
    See "systemctl status nvidia-persistenced.service" and "journalctl -xeu nvidia-persistenced.service" for details.

    No containers need to be restarted.

    No user sessions are running outdated binaries.

    No VM guests are running outdated hypervisor (qemu) binaries on this host.
   ```
   Still works, no issues observed.
8. Install Pytorch:
9. Install tmux:

### CPU-only
### GPU
