import os
from pathlib import Path
import pandas as pd
from ase.data import atomic_masses, atomic_numbers
import yaml
from .jobutils import *
from statebatch.build import build
import pprint

class Batch:
    def __init__(self, yaml_f, debug=False):
        with open (yaml_f, 'r') as f:
            config = yaml.safe_load(f)

        self.system_spec = config.get('system_spec')
        self.comp_spec   = config.get('comp_spec')
        self.dft_spec    = config.get('dft_spec')

        df          = pd.read_csv(self.comp_spec.get('csv_loc'))
        self.atoms_to_run = df.to_dict(orient='index')

        # Initializing class objects
        self.jobinfo = {}

        # Set the debug mode
        self.debug = debug
    
    def prerun(self, make_jobscript=None):
        def manage_system_params():
            for param in self.system_spec.get('fix_params'):
                for idx in range(len(self.atoms_to_run)):
                    self.atoms_to_run[idx][param] = self.system_spec.get('fix_params').get(param)
            
        def get_dft_params(atom_to_run, dft_spec, comp_spec):
            input_data = dft_spec.get('fix_params').copy()
            for param in dft_spec.get('vary_params'):
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
                os.symlink(os.path.join(comp_spec.get('pseudo_loc'), pseudo_file), pseudolink)
                pseudo_input = [element, atomic_masses[atomic_numbers[element]],pseudo_file]
                pseudo_list.append(pseudo_input)
            input_data['PSEUDOS'] = pseudo_list
            if input_data["XCTYPE"] in ['vdw-df', 'vdw-df2', 'rev-vdw-df2', 'optb86b-vdw']:
                input_data["VDW-DF"] =  {"QCUT": 10, "NQ": "20"}
            return input_data
        
        def link():
            pwlink = Path(os.path.join(os.getcwd(), self.comp_spec.get('pw_name')))
            pwlink.unlink(missing_ok=True)
            if os.path.exists(pwlink):
                os.remove(pwlink)
            os.symlink(os.path.join(self.comp_spec.get('pw_loc'), self.comp_spec.get('pw_name')), pwlink)

        manage_system_params()
        for idx, atom_to_run in enumerate(self.atoms_to_run.values(), start=1):
            cwd = os.getcwd()
            dirname = self.comp_spec.get('prefix')+str('{:04d}'.format(idx))
            os.makedirs(dirname, exist_ok=True)
            os.chdir(dirname)
            link()

            input_data = get_dft_params(atom_to_run, self.dft_spec, self.comp_spec)
            if self.debug:
                pprint.pprint(input_data)

            _, input_file, output_file = build(atom_to_run=atom_to_run,
                                               input_data=input_data,
                                               system_spec=self.system_spec,
                                               calc_name=self.dft_spec.get("name"),
                                               file_prefix=self.comp_spec.get("file_prefix")
                                               )
            os.chdir('../')

            # Save jobinfo
            self.jobinfo[idx] = {'idx':idx,
                                 'cwd':os.path.join(cwd,dirname),
                                 'input_file':input_file,
                                 'output_file': output_file}

        if make_jobscript is not None:
            write_jobscript(batch_obj=self, jobinfo=self.jobinfo, jobopt=make_jobscript)
