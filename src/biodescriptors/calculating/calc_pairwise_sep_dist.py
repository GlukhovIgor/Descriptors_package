import itertools
import numpy as np
import pandas as pd

from biodescriptors.calculating.calc_COM_helix import _calc_COM_helix
from biodescriptors.calculating import utils


def _calc_pairwise_sep_dist(chain, ref):
    """Calculate separation distance between every helix"""
    
    # Calculate centers of mass of every helix
    hel_COMs = list(itertools.chain(*_calc_COM_helix(chain, ref)))
    # Calculate pairwise separations between every COM's of every helix as vector distance
    pairwise_seps = []
    for j in range(0, len(hel_COMs) - 1):
        result = []
        for i in range(j + 1, len(hel_COMs)):
            vect = np.array(hel_COMs[i]) - np.array(hel_COMs[j])
            result.append(np.sqrt(np.sum(vect ** 2)))
        pairwise_seps.append(result)
    return pairwise_seps


def calc_pairwise_sep_dist(pdb_file, ref):
    """Calculate separation distance between every helix"""
    _, _, _, chain, _ = utils.get_model_and_structure(pdb_file)
    return _calc_pairwise_sep_dist(chain, ref)


def pairwise_sep_dist_to_pandas(pdb_file, ref, protein_name=None):
    """Putting separation distance between every helix in pandas dataframe."""
    cols_pairwise = ['prot_name'] + ['PairwiseSep H' + str(i) + '-H' + str(j) for i in range(1, 14) for j in range(i+1, 14)]
    df_pairseps = pd.DataFrame(columns=cols_pairwise)
    pairseps = None
    try:
        pairseps = calc_pairwise_sep_dist(pdb_file, ref)
    except:
        KeyError
        print('KeyError while calculating pairsep')
    data_pairseps = [protein_name]
    if pairseps is not None:
        for elem in pairseps:
            for dist in elem:
                data_pairseps.append(dist)
    df_pairseps = df_pairseps.append(pd.Series(data_pairseps, index=cols_pairwise[0:len(data_pairseps)]), ignore_index=True)
    return df_pairseps
