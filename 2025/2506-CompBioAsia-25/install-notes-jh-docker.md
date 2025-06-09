# Installation for Amber MD Prep

See https://github.com/CompBioAsia/Amber-md-prep (private, requires access)

## Required Tools

For `cba_tools.py`:
- [mdtraj](https://mdtraj.org/1.9.3/installation.html)
- [pdbfixer](https://github.com/openmm/pdbfixer) -- might want to use the following installation, as standard pdbfixer uses conda: https://pypi.org/project/black-pdbfixer/ (installs 1.9 instead of latest (1.11))
- [openmm](http://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm) -- NOTE: install GPU-ready version
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