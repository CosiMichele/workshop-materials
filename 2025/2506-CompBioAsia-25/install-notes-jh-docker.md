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
See `Dockerfile` for more details. Docker container was built with `docker build .` and pushed to Dockerhub as `cosimichele/jupyter-amber:gpu-250608`, where it can be pulled and deployed with [CACAO](https://cacao.jetstream-cloud.org/).

## Running Training

Notes:
- 5th box is `add_h('abl_imatinib_heavy.pdb', 'abl_imatinib_amber.pdb', chimerax='chimerax', mode='amber')`. Should be `add_h('abl_imatinib_heavy.pdb', 'abl_imatinib_amber.pdb', chimera='chimerax', mode='amber')`
- `add_h` function in `cba_tools.py` is buggy. Suggested workaround for `add_h`:
    ```
    def add_h(inpdb, outpdb, chimera='chimerax', mode='amber'):
        '''
        Add hydrogen atoms to a PDB format file, using ChimeraX (not Chimera)

        Args:
            inpdb (str): name of input PDB file
            outpdb (str): name of output PDB file
            chimera (str): command to invoke ChimeraX
            mode (str): Adjust residue names for chosen software
                        (only 'amber' is currently supported).
        '''
        check_available(chimera)
        fh = FileHandler()
        script = fh.create('script')

        # ChimeraX-compatible syntax
        script.write_text(
            'open infile.pdb\n'
            'addh\n'
            'save outfile.pdb format pdb\n'
            'exit\n'
        )

        addh = SubprocessKernel(f"{chimera} --nogui < script")
        addh.set_inputs(['script', 'infile.pdb'])
        addh.set_outputs(['outfile.pdb'])
        addh.set_constant('script', script)

        infile = fh.load(inpdb)
        addh_result = addh.run(infile)

        if addh_result is None:
            raise RuntimeError("ChimeraX failed to generate outfile.pdb")

        if mode == 'amber':
            check_available('pdb4amber')
            pdb4amber = SubprocessKernel('pdb4amber -i infile.pdb -o outfile.pdb')
            pdb4amber.set_inputs(['infile.pdb'])
            pdb4amber.set_outputs(['outfile.pdb'])
            amberpdb = pdb4amber.run(addh_result)
            amberpdb.save(outpdb)
        else:
            addh_result.save(outpdb)

        print(f'Hydrated structure written to {outpdb}')
    ```
    - Code runs. Need to check output