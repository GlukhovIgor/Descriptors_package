import numpy as np
import pandas as pd

from biodescriptors.calculating import constraints
from biodescriptors.calculating import utils


def _calc_COM_helix(chain, ref):
    """Calculate center of mass for every helix in PDB structure"""

    helices = ref
    hel_COM = []

    for elem in helices:
        helix_mass = 0
        weighted_coord = list()
        helix_com = list()

        for res in elem:
            residue = chain[res]

            for atom in residue.get_atoms():
                weight = constraints.ATOMIC_WEIGHTS[atom.get_name()[0]]
                helix_mass += weight
                weighted_coord.append([coord * weight for coord in list(atom.get_coord())])

        helix_com.append([coord / helix_mass for coord in np.sum(weighted_coord, axis=0)])
        hel_COM.append(helix_com)

    return hel_COM


def calc_COM_helix(pdb_file, ref):
    """Calculate center of mass for every helix in PDB structure"""
    # Initialize PDB structure
    _, _, _, chain, _ = utils.get_model_and_structure(pdb_file)
    return _calc_COM_helix(chain, ref)


def COM_helix_to_pandas(pdb_file, ref, protein_name=None):
    """Putting center of mass for every helix in PDB structure in pandas dataframe."""
# Not implemented yet
    return None
