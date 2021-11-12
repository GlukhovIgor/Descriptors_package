import numpy as np
import pandas as pd

from biodescriptors.calculating import constraints
from biodescriptors.calculating import utils


def _calc_COM_protein(atom_struct):
    """Calculate protein's center of mass"""

    # Calculate product of atom's coordintates and atomic weights of each atom
    atoms = [([coord * constraints.ATOMIC_WEIGHTS[atom.get_name()[0]] for coord in list(atom.get_coord())]) for atom in atom_struct]

    # Calculate total mass of PDB structure
    total_mass = sum([constraints.ATOMIC_WEIGHTS[atom.get_name()[0]] for atom in atom_struct])

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
