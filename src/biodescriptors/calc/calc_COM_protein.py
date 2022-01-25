import numpy as np
import pandas as pd

from biodescriptors.calc import constraints
from biodescriptors.calc import utils


def _calc_COM_protein(atom_struct):
    """Calculate protein's center of mass""" 
    atoms = []
    total_mass = 0
    for atom in atom_struct:
        
        # Calculate product of atom's coordintates and atomic weights of each atom
        single_atom_com = []
        for coord in list(atom.get_coord()):
            single_atom_com.append(coord * constraints.ATOMIC_WEIGHTS[atom.get_name()[0]])
        atoms.append(single_atom_com)
        
        # Calculate total mass of PDB structure
        total_mass += constraints.ATOMIC_WEIGHTS[atom.get_name()[0]]

    # Calculate an protein's center of mass
    COM = [coord / total_mass for coord in np.sum(atoms, axis=0)]
    return COM


def calc_COM_protein(pdb_file):
    """Calculate protein's center of mass"""
    # Initialize PDB structure
    _, _, _, _, atom_struct = utils.get_model_and_structure(pdb_file)
    return _calc_COM_protein(atom_struct)


def COM_protein_to_pandas(pdb_file, protein_name=None):
    """Putting protein's center of mass in pandas dataframe."""
# Not implemented
    return None
