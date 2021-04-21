from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from pymol import cmd
from parameters import Parameters

def remove_chains(structure, chain_to_keep):


    model = structure[0]

    chains_to_remove = []
    for chain in model.child_list:
        if chain.id != chain_to_keep:
            chains_to_remove.append(chain.id)

    for id in chains_to_remove:
        model.detach_child(id)

    return structure



def remove_extra_atoms(structure, parameters):
    model = structure[0]
    chain = model.child_list[0]
    residues_to_remove = []
    for residue in chain.child_list:
        if not parameters.has_residue(residue.resname):
            residues_to_remove.append(residue.id)
            continue
        res = residue.resname
        atoms_to_remove = []
        for atom in residue.child_list:
            atom_name = atom.name
            if atom_name == "OP1":
                atom_name = "O1P"
            if atom_name == "OP2":
                atom_name = "O2P"
            if not parameters.has_atom(res, atom_name):
                atoms_to_remove.append(atom.id)

        for atom in atoms_to_remove:
            residue.detach_child(atom)

        residue.sort()
    for residue in residues_to_remove:
        chain.detach_child(residue)
    return structure


def clean_pdb(file_input, file_output, chain_to_keep, parameters):
    parser = MMCIFParser()

    structure = parser.get_structure(file_input[:-4].upper(), file_input)

    structure = remove_chains(structure, chain_to_keep)
    structure = remove_extra_atoms(structure, parameters)

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(file_output)
