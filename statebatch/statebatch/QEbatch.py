import os, sys
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
from ase.calculators.espresso import Espresso
from ase import Atoms
from ase.build import molecule, bulk, fcc100, fcc110, fcc111, bcc100, bcc110, bcc111, add_adsorbate
from ase.data import atomic_masses, atomic_numbers
from ase.db import connect
import yaml
from pprint import pprint
from .jobutils import *

class Batch:
    """Batch wrapper

    Parameters
    ----------
    yaml_f : str, path object or file-like object
        YAML file path
    """
    def __init__(self, yaml_f):
        # Read YAML file
        with open (yaml_f, 'r') as f:
            config = yaml.safe_load(f)

        # Read YAML input general components
        self.system_spec = config.get('system_spec')
        self.comp_spec   = config.get('comp_spec')
        self.dft_spec    = config.get('dft_spec')

        # Read CSV file (system and parameter list)
        df = pd.read_csv(self.comp_spec.get('csv_loc'))
        self.atoms_to_run = df.to_dict(orient='index')

        # Initializing class objects
        self.jobinfo = {}

    def prerun(self, make_jobscript=None):
        """Prepares work directory and run files"""
        def manage_system_params():
            """Passes system `fix_params` to atoms_to_run dictionary"""
            for param in self.system_spec.get('fix_params'):
                for idx in range(len(self.atoms_to_run)):
                    self.atoms_to_run[idx][param] = self.system_spec.get('fix_params').get(param)

        def get_dft_params(atom_to_run):
            """Distributes dft parameters to calculator input"""
            # Initialize input_data
            input_data = self.dft_spec.get('fix_params').copy()

            # Replace template params with vary_params
            for param in self.dft_spec.get('vary_params'):
                input_data[param] = atom_to_run[param]

            # Pseudopotential preparation
            # -- Atomic mass not implemented here --
            pseudos = atom_to_run['PSEUDOS']
            pseudopotentials = {}
            for pseudo_element_pair in pseudos.split('|'):
                element = pseudo_element_pair.split('@')[0]
                pseudo_file = pseudo_element_pair.split('@')[1]
                pseudopotentials[element] = pseudo_file
            input_data.pop('PSEUDOS', 'not_found')

            input_data['pseudopotentials'] = pseudopotentials # **
            input_data['PSEUDO_DIR'] = self.comp_spec.get('pseudo_loc')

            # XC functional preparation
            if {'input_dft', 'vdw_corr'} <=  input_data.keys():
                raise Exception('Use either input_dft or vdw_corr (not together)')
            if input_data['XCTYPE'] in ['grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d',
                                        'grimme-d3', 'Grimme-D3', 'DFT-D3', 'dft-d3',
                                        'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler',
                                        'MBD', 'mbd', 'many-body-dispersion', 'mbd_vdw',
                                        'XDM', 'xdm']:
                input_data['VDW_CORR'] = input_data['XCTYPE']
            else:
                input_data['INPUT_DFT'] = input_data['XCTYPE']
            input_data.pop('XCTYPE', 'not_found')

            # K-points preparation
            kpts = input_data['KPTS'].copy()
            input_data['kpts'] = kpts

            return (input_data)

        def build(atom_to_run):
            """Build

            Atomic structure builder and input file writer

            Parameters
            ----------
            atom_to_run : dict
                Atomic structure dictionary
            """
            if (self.system_spec.get('type') == 'Atom'):
                _atoms       = atom_to_run['Species']
                atoms_obj = Atoms(_atoms)
                atoms_obj.set_cell(atom_to_run['Vacuum']*np.identity(3))
            elif (self.system_spec.get('type') == 'Molecule'):
                _atoms = atom_to_run['Molecule']
                atoms_obj = molecule(_atoms)
                atoms_obj.set_cell(atom_to_run['Vacuum']*np.identity(3))
            elif (self.system_spec.get('type') == 'Bulk'):
                _atoms = atom_to_run['Bulk']
                crystalstructure = atom_to_run['Crystalstructure']
                atoms_obj = bulk(_atoms, crystalstructure=crystalstructure)
            elif (self.system_spec.get('type') == 'Surface' or self.system_spec.get('type') == 'Adsorption'):
                _atoms = atom_to_run['Surface']
                crystalstructure = atom_to_run['Crystalstructure']
                facet = str(atom_to_run['Facet'])
                size = eval(atom_to_run['Size'])
                vacuum = 0.5*atom_to_run['Vacuum']
                if (crystalstructure+facet == 'fcc100'):
                    atoms_obj = fcc100(_atoms, size=size, vacuum=vacuum)
                elif (crystalstructure+facet == 'fcc110'):
                    atoms_obj = fcc110(_atoms, size=size, vacuum=vacuum)
                elif (crystalstructure+facet == 'fcc111'):
                    atoms_obj = fcc111(_atoms, size=size, vacuum=vacuum)
                elif (crystalstructure+facet == 'bcc100'):
                    atoms_obj = bcc100(_atoms, size=size, vacuum=system_vacuum)
                elif (crystalstructure+facet == 'bcc110'):
                    atoms_obj = bcc110(_atoms, size=size, vacuum=system_vacuum)
                elif (crystalstructure+facet == 'bcc111'):
                    atoms_obj = bcc111(_atoms, size=size, vacuum=system_vacuum)

                if (self.system_spec.get('type') == 'Adsorption'):
                    _adsorbate = atom_to_run['Adsorbate']
                    adsorbate_obj = molecule(_adsorbate)
                    add_adsorbate(atoms_obj, adsorbate_obj, height = atom_to_run['Height'], position = atom_to_run['Site'])

            # Finalize input_data and input file
            input_data = get_dft_params(atom_to_run)
            if self.comp_spec.get('fileprefix'):
                label = self.comp_spec.get('fileprefix')
            else:
                label = f"state"

            input_file, output_file = f"{label}.in", f"{label}.out"
            atoms_obj.calc = STATE(label=label, input_data=input_data)
            atoms_obj.calc.write_input(atoms_obj)
            atoms_obj.write(f"{_atoms}.xyz")        #For using prefix: atoms_obj.write(f"{label}.xyz")

            return (atoms_obj, input_file, output_file)

        # Run prerun for all systems in CSV
        manage_system_params()
        for idx in range(len(self.atoms_to_run)):

            cwd = os.getcwd()
            dirname = self.comp_spec.get('prefix')+str('{:04d}'.format(idx))
            os.makedirs(dirname, exist_ok=True)
            os.chdir(dirname)
            _, input_file, output_file = build(self.atoms_to_run[idx])
            os.chdir('../')

            # Save jobinfo
            self.jobinfo[idx] = {'idx':idx,
                                 'cwd':os.path.join(cwd,dirname),
                                 'input_file':input_file,
                                 'output_file': output_file}

        if make_jobscript is not None:
            write_jobscript(batch_obj=self, jobinfo=self.jobinfo, jobopt=make_jobscript)
