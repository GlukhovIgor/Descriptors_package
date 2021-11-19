import numpy as np
import pandas as pd

from biodescriptors.calculating.calc_COM_protein import _calc_COM_protein
from biodescriptors.calculating import utils


def _calc_COM_Calpha_angles(chain, atom_struct, ref):
    """Calculate angles between protein's center of mass and alpha carbon atom of every helix"""
    
    # Initialize PDB structure
    protCOM = _calc_COM_protein(atom_struct)
    helix_content = ref 

    # Calculate coordinates of alpha carbon atom of every helix
    CA_coords = []
    for elem in helix_content:
        helices = [list(elem)]
        for elem in helices:
            coord = []
            for res in elem:
                residue = chain[res]
                for atom in residue.get_atoms():
                    if atom.get_name() == 'CA':
                        coord.append(atom.get_coord())
            CA_coords.append(coord)

    # Calculate angles between protein's center of mass and helices of every alpha carbon
    angles = []
    for elem in CA_coords:
        # Calculate it as angle formed by two vectors
        OHel1 = elem[0] - np.array(protCOM)
        OHel2 = elem[1] - np.array(protCOM)

        OHels = np.dot(OHel1,OHel2)

        OHel1abs = np.linalg.norm(OHel1)
        OHel2abs = np.linalg.norm(OHel2)
        angles.append(np.degrees(np.arccos((OHels / (OHel1abs * OHel2abs)))))
    return angles


def calc_COM_Calpha_angles(pdb_file, ref):
  #rename
    """Calculate angles between protein's center of mass and alpha carbon atom of every helix"""
    _, _, _, chain, atom_struct = utils.get_model_and_structure(pdb_file)
    return _calc_COM_Calpha_angles(chain, atom_struct, ref)


def COM_Calpha_angles_to_pandas(pdb_file, ref, protein_name=None):
    """Putting angles between protein's center of mass and alpha carbon atom of every helix in pandas dataframe."""
    cols_angle = ['prot_name'] + ['AngleCOM H' + str(elem) for elem in range(1, 14)]
    df_alphaangle = pd.DataFrame(columns=cols_angle)
    alpha_angle = None
    try:
        alpha_angle = calc_COM_Calpha_angles(pdb_file, ref)
    except KeyError:
        print('KeyError while calculating alpha angle')
    data_alphaagnle = [protein_name]
    if alpha_angle is not None:
        for elem in alpha_angle:
            data_alphaagnle.append(elem)
    df_alphaangle = df_alphaangle.append(pd.Series(data_alphaagnle, index=cols_angle[0:len(data_alphaagnle)]), ignore_index=True)
    return df_alphaangle
