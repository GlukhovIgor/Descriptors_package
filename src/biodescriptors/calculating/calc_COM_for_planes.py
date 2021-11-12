import numpy as np
import pandas as pd

from biodescriptors.calculating import constraints
from biodescriptors.calculating import utils


def _calc_COM_for_planes(chain, helices):
    """Calculate center of mass for every "sandwich layer" of VDR structure. Helices - list of layer's helices"""
    
    # Calculate center of mass for every helix
    hel_COM = []
    
    for elem in helices:
        helix_content = [list(elem)]
        helix_mass = 0
        weighted_coord = list()
        helix_com = list()

        for elem in helix_content:
            for res in elem:
                residue = chain[res]

                for atom in residue.get_atoms():
                    weight = constraints.ATOMIC_WEIGHTS[atom.get_name()[0]]
                    helix_mass += weight # Calculate helix total mass
                    # Calculate product of atom coordinate and weight
                    weighted_coord.append([coord * weight for coord in list(atom.get_coord())])
            # Calculate helix center of mass
            helix_com.append([coord / helix_mass for coord in np.sum(weighted_coord, axis=0)])
            hel_COM.append(helix_com)
    return hel_COM


def calc_COM_for_planes(pdb_file, helices):
    """Calculate center of mass for every "sandwich layer" of VDR structure. Helices - list of layer's helices."""

    _, _, _, chain, _ = utils.get_model_and_structure(pdb_file)
    return _calc_COM_for_planes(chain, helices)


def COM_for_planes_to_pandas(pdb_file, helices, protein_name=None):
    """Putting center of mass for every "sandwich layer" of VDR structure in pandas dataframe."""
# not implemented yet
    return None
