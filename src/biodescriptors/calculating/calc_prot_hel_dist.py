import numpy as np
import pandas as pd

from biodescriptors.calculating.calc_COM_helix import _calc_COM_helix
from biodescriptors.calculating.calc_COM_protein import _calc_COM_protein
from biodescriptors.calculating import utils


def _calc_prot_hel_dist(chain, atom_struct, ref):
    """Calculate distance between protein's center of mass and between every helix's center of mass"""
    
    # Calculate protein's center of mass and centers of helices masses
    hel_COMs = _calc_COM_helix(chain, ref)
    prot_COM = _calc_COM_protein(atom_struct)
    # Calculate distance between protein's COM and helices COMs as vector distance
    prot_hel_dists = []
    for i in range(0, len(hel_COMs)):
        vect = np.array(hel_COMs[i]) - np.array(prot_COM)
        prot_hel_dists.append(np.sqrt(np.sum(vect ** 2)))
    return prot_hel_dists


def calc_prot_hel_dist(pdb_file, ref):
    """Calculate distance between protein's center of mass and between every helix's center of mass"""
    _, _, _, chain, atom_struct = utils.get_model_and_structure(pdb_file)
    return _calc_prot_hel_dist(chain, atom_struct, ref)


def prot_hel_dist_to_pandas(pdb_file, ref, protein_name=None):
    """Putting distance between protein's center of mass and between every helix's center of mass in pandas dataframe."""
    cols_protheldist = ['prot_name'] + ['Dist prot-H' + str(elem) for elem in range(1, 14)]
    df_prothel = pd.DataFrame(columns=cols_protheldist)
    prothel = None
    try:
        prothel = calc_prot_hel_dist(pdb_file, ref)
    except:
        KeyError
        print('KeyError while calculating prot-helix distance')
    data_prothel = [protein_name]
    if prothel is not None:
        for elem in prothel:
            data_prothel.append(elem)
    df_prothel = df_prothel.append(pd.Series(data_prothel, index=cols_protheldist[0:len(data_prothel)]), ignore_index=True)
    return df_prothel
