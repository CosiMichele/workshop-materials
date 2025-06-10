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

All VMs will have the following:
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) & [Mamba](https://mamba.readthedocs.io/en/latest/index.html).
- Python (3.11)
- Pandas
- Matplotlib
- Numpy
- Gemmi
- Pytorch (note: requires GPU VM for building)

For CPU only VMs, these sofware will be installed on top of universal software:
- HMMER
- MMSeq2
- nail

For GPU VMs, these sofware will be installed on top of universal software:
- AMBERTools25
- AMBER24
- AutoDock Vina

TBD:
- Alphafold2 
- Alphafold3

