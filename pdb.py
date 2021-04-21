from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import PDBIO

from Bio.PDB import Model
import numpy as np
# This file handles everything related to the pdb structure such as atom indices, keeping track of position and velocity vectors etc

class PDB(object):



    def __init__(self, file):
        parser = MMCIFParser()
        self.structure = parser.get_structure(file[:4].upper(), file)
        self.initial_conditions()

    # create the position + velocity vector from the structure
    def initial_conditions(self):
        model = self.structure[0]
        chain = model.child_list[0]

        #positions = np.empty(( , 3))
        atoms = []
        coords = []

        for residue in chain.get_residues():
            for atom in residue.get_atoms():
                atoms.append(atom)
                coords.append(atom.get_coord())

        self.positions = np.asarray(coords)
        self.atoms = atoms

    # sets new positions and add new model
    def set_new_positions(self, positions):
        self.positions = positions
        new_model = Model.Model(len(self.structure.child_list), len(self.structure.child_list) + 1)
        chain = self.structure[0].child_list[0].copy()
        new_model.add(chain)

        self.structure.add(new_model)
        chain = new_model.child_list[0]
        counter = 0
        for residue in chain.get_residues():
            for atom in residue.get_atoms():
                atom.set_coord(positions[counter])
                counter += 1


    def save(self, fileout):

        #io = PDBIO()
        #io.set_structure(self.structure)
        #io.save(fileout)
        io = MMCIFIO()
        io.set_structure(self.structure)
        io.save(fileout)

# how to store atom info?
# need to check if atoms in same residue
# only calculate 1 half for potentials, symmetric

