from statebatch.build import build_atoms_structure


def test_build_atoms_structure():
    system_spec = {
        "type": "Atom", 
        "fix_params": {"Vacuum": 30.0, "Alat": "from_ase", "Size": "(1,1,1)"}
    }

    atom_to_run = {
        "Species": "H",
        "PSEUDOS": "H@h_pbe_v1.4.uspp.F.UPF",
        "XCTYPE": "pbe",
    }

    for param in system_spec.get('fix_params'):
        atom_to_run[param] = system_spec.get('fix_params').get(param)

    atoms_obj, _atoms = build_atoms_structure(atom_to_run, system_spec)

    # Check if an atoms object has been created
    assert len(atoms_obj) == 1
