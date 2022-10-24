from statebatch.build import build_atoms_structure


def test_build_atoms_structure():
    system_spec = {
        "type": "Adsorption",
        "fix_params": {
            "Vacuum": 30.0,
            "Alat": "from_ase",
            "Size": "(3,3,3)",
            "Crystalstructure": "fcc",
            "Facet": "111",
            "Height": 1.8
        }
    }

    atom_to_run = {
        "Surface": "Cu",
        "Adsorbate": "O",
        "PSEUDOS": "Cu@pot.Cu_pbe1|C@pot.C_pbe1|O@pot.O_pbe1",
        "XCTYPE": "ggapbe",
        "Site": "ontop"
    }

    for param in system_spec.get('fix_params'):
        atom_to_run[param] = system_spec.get('fix_params').get(param)

    atoms_obj, _atoms = build_atoms_structure(atom_to_run, system_spec)

    # Check if an atoms object has been created
    assert len(atoms_obj) == 28
