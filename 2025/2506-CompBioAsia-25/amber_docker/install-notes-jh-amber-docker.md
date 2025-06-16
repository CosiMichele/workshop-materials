# Installation for Amber MD Prep

See https://github.com/CompBioAsia/Amber-md-prep (private, requires access)

## Required Tools

For `cba_tools.py`:
- [mdtraj](https://mdtraj.org/1.9.3/installation.html)
- [pdbfixer](https://github.com/openmm/pdbfixer) -- might want to use the following installation, as standard pdbfixer uses conda: https://pypi.org/project/black-pdbfixer/ (installs 1.9 instead of latest (1.11))
- [openmm](http://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm)
- [crossflow](https://pypi.org/project/crossflow/)
- [numpy](https://pypi.org/project/numpy/)
- functools (installed in standard lib)
- shutil (installed in standard lib)

For `protein_ligand_prep.ipynb`:
- mdtraj (link above)
- [nglview](https://github.com/nglviewer/nglview)
- cba_tools (`cba_tools.py`?)
- openmm (link above)

Required for functionality:
- conda (mamba)
    - cudatoolkit
    - ipywidgets=7.7.2 
    - gitpython 
    - nbconvert 
    - nbformat 
    - requests 
    - packaging 
    - rich 
    - jax 
    - parmed 
    - netCDF4 
    - openbabel 
    - widgetsnbextension 
    - jupyter_contrib_nbextensions 
    - nodejs
    - ambertools
    - rdkit
    - biobb_amber (bioconda)

Additional functionality:
- `entry.sh`
- `jupyter_notebook_config.json`

## Building Docker
See `Dockerfile` for more details. Docker container was built with `docker build .` and pushed to Dockerhub as **`cosimichele/jupyter-amber:gpu-250609`**, where it can be pulled and deployed with [CACAO](https://cacao.jetstream-cloud.org/).

## Running Training

Notes:
- (2025-06-15) Running **protein-ligand-prep**,  cell 3 looks for file `ligand_PH7.pdb`, should be `ligand_pH7.pdb`
- (2025-06-15) **openMM-md-run** runs end to end.