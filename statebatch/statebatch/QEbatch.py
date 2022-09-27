import os, sys
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
from state_interface.state import STATE
from ase import Atoms
from ase.build import molecule, bulk, fcc100, fcc110, fcc111, bcc100, bcc110, bcc111, add_adsorbate
from ase.data import atomic_masses, atomic_numbers
from ase.db import connect
import yaml

class Batch:
    def __init__(self, yaml_f):
        with open (yaml_f, 'r') as f:
            config = yaml.safe_load(f)

        self.system_spec = config.get('system_spec')
        self.comp_spec   = config.get('comp_spec')
        self.dft_spec    = config.get('dft_spec')
        
        df          = pd.read_csv(self.comp_spec.get('csv_loc'))
        self.atoms_to_run = df.to_dict(orient='index')

    def prerun(self):
        def manage_system_params():
            for param in self.system_spec.get('fix_params'):
                for idx in range(len(self.atoms_to_run)):
                    self.atoms_to_run[idx][param] = self.system_spec.get('fix_params').get(param)
        def get_dft_params(atom_to_run):
            input_data = self.dft_spec.get('fix_params').copy()
            for param in self.dft_spec.get('vary_params'):
                input_data[param] = atom_to_run[param]
            pseudos = atom_to_run['PSEUDOS']
            pseudo_list = []
            for pseudo_element_pair in pseudos.split('|'):
                element = pseudo_element_pair.split('@')[0]
                pseudo_file = pseudo_element_pair.split('@')[1]
                pseudolink = Path(os.path.join(os.getcwd(), pseudo_file))
                pseudolink.unlink(missing_ok=True)
                if os.path.exists(pseudolink):
                    os.remove(pseudolink)
                os.symlink(os.path.join(self.comp_spec.get('pseudo_loc'), pseudo_file), pseudolink)
                pseudo_input = [element, atomic_masses[atomic_numbers[element]],pseudo_file]
                pseudo_list.append(pseudo_input)
            input_data['PSEUDOS'] = pseudo_list
            if input_data["XCTYPE"] in ['vdw-df', 'vdw-df2', 'rev-vdw-df2', 'optb86b-vdw']:
                input_data["VDW-DF"] =  {"QCUT": 10, "NQ": "20"}
            return (input_data)
        def link():
            pwlink = Path(os.path.join(os.getcwd(), self.comp_spec.get('pw_name')))
            pwlink.unlink(missing_ok=True)
            if os.path.exists(pwlink):
                os.remove(pwlink)
            os.symlink(os.path.join(self.comp_spec.get('pw_loc'), self.comp_spec.get('pw_name')), pwlink)
        def build(atom_to_run):
            if (self.system_spec.get('type') == 'Atom'):
                _atoms       = atom_to_run['Species']
                self.atoms_obj = Atoms(_atoms)
                self.atoms_obj.set_cell(atom_to_run['Vacuum']*np.identity(3))
            elif (self.system_spec.get('type') == 'Molecule'):
                _atoms = atom_to_run['Molecule']
                self.atoms_obj = molecule(_atoms)
                self.atoms_obj.set_cell(atom_to_run['Vacuum']*np.identity(3))
            elif (self.system_spec.get('type') == 'Bulk'):
                _atoms = atom_to_run['Bulk']
                crystalstructure = atom_to_run['Crystalstructure']
                self.atoms_obj = bulk(_atoms, crystalstructure=crystalstructure)
            elif (self.system_spec.get('type') == 'Surface' or self.system_spec.get('type') == 'Adsorption'):
                _atoms = atom_to_run['Surface']
                crystalstructure = atom_to_run['Crystalstructure']
                facet = str(atom_to_run['Facet'])
                size = eval(atom_to_run['Size'])
                vacuum = 0.5*atom_to_run['Vacuum']
                if (crystalstructure+facet == 'fcc100'):
                    self.atoms_obj = fcc100(_atoms, size=size, vacuum=vacuum)
                elif (crystalstructure+facet == 'fcc110'):
                    self.atoms_obj = fcc110(_atoms, size=size, vacuum=vacuum)
                elif (crystalstructure+facet == 'fcc111'):
                    self.atoms_obj = fcc111(_atoms, size=size, vacuum=vacuum)
                elif (crystalstructure+facet == 'bcc100'):
                    self.atoms_obj = bcc100(_atoms, size=size, vacuum=system_vacuum)
                elif (crystalstructure+facet == 'bcc110'):
                    self.atoms_obj = bcc110(_atoms, size=size, vacuum=system_vacuum)
                elif (crystalstructure+facet == 'bcc111'):
                    self.atoms_obj = bcc111(_atoms, size=size, vacuum=system_vacuum)

                if (self.system_spec.get('type') == 'Adsorption'):
                    _adsorbate = atom_to_run['Adsorbate']
                    self.adsorbate_obj = molecule(_adsorbate)
                    add_adsorbate(self.atoms_obj, self.adsorbate_obj, height = atom_to_run['Height'], position = atom_to_run['Site'])


            input_data = get_dft_params(atom_to_run)
            label = f"{_atoms}"
            input_file, output_file = f"{label}.in", f"{label}.out"
            self.atoms_obj.calc = STATE(label=label, input_data=input_data)
            self.atoms_obj.calc.write_input(self.atoms_obj)
            self.atoms_obj.write(f"{label}.xyz")
        
        manage_system_params()
        for idx in range(len(self.atoms_to_run)):
            dirname = self.comp_spec.get('prefix')+str(idx)
            os.makedirs(dirname, exist_ok=True)
            os.chdir(dirname)
            link()
            build(self.atoms_to_run[idx])
            os.chdir('../')