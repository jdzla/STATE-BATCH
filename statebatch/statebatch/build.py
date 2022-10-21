from ase import Atoms
import numpy as np
from ase.build import (
    molecule,
    bulk,
    fcc100,
    fcc110,
    fcc111,
    bcc100,
    bcc110,
    bcc111,
    add_adsorbate,
)


def build_atoms_structure(atom_to_run, system_spec):
    """Build

    Atomic structure builder and input file writer

    Parameters
    ----------
    atom_to_run : dict
        Atomic structure dictionary
    system_spec :

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
        if crystalstructure + facet == "fcc100":
            atoms_obj = fcc100(_atoms, size=size, vacuum=vacuum)
        elif crystalstructure + facet == "fcc110":
            atoms_obj = fcc110(_atoms, size=size, vacuum=vacuum)
        elif crystalstructure + facet == "fcc111":
            atoms_obj = fcc111(_atoms, size=size, vacuum=vacuum)
        elif crystalstructure + facet == "bcc100":
            atoms_obj = bcc100(_atoms, size=size, vacuum=vacuum)
        elif crystalstructure + facet == "bcc110":
            atoms_obj = bcc110(_atoms, size=size, vacuum=vacuum)
        elif crystalstructure + facet == "bcc111":
            atoms_obj = bcc111(_atoms, size=size, vacuum=vacuum)

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
