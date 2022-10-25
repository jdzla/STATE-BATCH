from ase import Atoms
import numpy as np
from ase.build import (
    molecule,
    bulk,
    add_adsorbate,
)
from state_interface.state import STATE
from ase.calculators.espresso import Espresso
import importlib


def build(atom_to_run, input_data, system_spec, calc_name, file_prefix=None):
    atoms_obj, _atoms = build_atoms_structure(atom_to_run, system_spec)
    if file_prefix:
        label = file_prefix
    else:
        if calc_name == "STATE":
            label = "state"
        if calc_name == 'QE':
            label = "espresso"

    if calc_name == "STATE":
        input_file, output_file = f"{label}.in", f"{label}.out"
        atoms_obj.calc = STATE(label=label, input_data=input_data)
    elif calc_name == "QE":
        input_file, output_file = f"{label}.pwi", f"{label}.pwo"
        atoms_obj.calc = Espresso(label=label, **input_data)
    atoms_obj.calc.write_input(atoms_obj)
    # For using prefix: atoms_obj.write(f"{label}.xyz")
    atoms_obj.write(f"{_atoms}.xyz")

    return (atoms_obj, input_file, output_file)


def build_atoms_structure(atom_to_run, system_spec):
    """Build

    Atomic structure builder

    Parameters
    ----------
    atom_to_run : dict
        Atomic structure dictionary
    system_spec : dict
        Information about the system to be calculated

    """
    if system_spec.get("type") == "Atom":
        _atoms = atom_to_run["Species"]
        atoms_obj = Atoms(_atoms)
        atoms_obj.set_cell(atom_to_run["Vacuum"] * np.identity(3))
    elif system_spec.get("type") == "Molecule":
        _atoms = atom_to_run["Molecule"]
        atoms_obj = molecule(_atoms)
        atoms_obj.set_cell(atom_to_run["Vacuum"] * np.identity(3))
    elif system_spec.get("type") == "Bulk":
        _atoms = atom_to_run["Bulk"]
        crystalstructure = atom_to_run["Crystalstructure"]
        atoms_obj = bulk(_atoms, crystalstructure=crystalstructure)
    elif (
        system_spec.get("type") == "Surface" or system_spec.get("type") == "Adsorption"
    ):
        _atoms = atom_to_run["Surface"]
        crystalstructure = atom_to_run["Crystalstructure"]
        facet = str(atom_to_run["Facet"])
        size = eval(atom_to_run["Size"])
        vacuum = 0.5 * atom_to_run["Vacuum"]
        surface_builder_function = getattr(
            importlib.import_module('ase.build'),
            crystalstructure + facet
        )
        atoms_obj = surface_builder_function(_atoms, size=size, vacuum=vacuum)

        if system_spec.get("type") == "Adsorption":
            _adsorbate = atom_to_run["Adsorbate"]
            adsorbate_obj = molecule(_adsorbate)
            add_adsorbate(
                atoms_obj,
                adsorbate_obj,
                height=atom_to_run["Height"],
                position=atom_to_run["Site"],
            )

    return atoms_obj, _atoms
