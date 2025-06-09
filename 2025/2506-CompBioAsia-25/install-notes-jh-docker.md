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
    - biobb_amber (bioconda)

Additional functionality:
- `entry.sh`
- `jupyter_notebook_config.json`

## Building Docker
See `Dockerfile` for more details. Docker container was built with `docker build .` and pushed to Dockerhub as `cosimichele/jupyter-amber:gpu-250608`, where it can be pulled and deployed with [CACAO](https://cacao.jetstream-cloud.org/).

## Running Training

Notes:
- 5th box is `add_h('abl_imatinib_heavy.pdb', 'abl_imatinib_amber.pdb', chimerax='chimerax', mode='amber')`. Should be `add_h('abl_imatinib_heavy.pdb', 'abl_imatinib_amber.pdb', chimera='chimerax', mode='amber')`
- `add_h` function in `cba_tools.py` is buggy. Suggested workaround?