import pandas as pd

from biodescriptors.calculating import utils


def _calc_len_of_hel(chain, ref):
    """Calculation of length of helices from structure"""
    
    helix_borders = ref

    # extracting vectors of coordinates for every border
    helices = {}
    for i in range(0, len(helix_borders)):
        helices[i] = [chain[res]['CA'].get_vector() for res in [helix_borders[i]][0]]

    for elem in helices:
        helices[elem] = (helices[elem][1]-helices[elem][0])

    lens_of_helices = []

    # calculation of length of helices using formula ((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)^0.5
    for el in helices:
        lens_of_helices.append((helices[el][0]**2 + helices[el][1]**2 + helices[el][2]**2) ** 0.5)

    return lens_of_helices


def calc_len_of_hel(pdb_file, ref):
    """Calculation of length of helices from structure"""
    _, _, _, chain, _ = utils.get_model_and_structure(pdb_file)
    return _calc_len_of_hel(chain, ref)


def len_of_hel_to_pandas(pdb_file, ref, protein_name=None):
    """Putting length of helices from structure in pandas dataframe."""
    cols_len = ['prot_name'] + [f'Length H{elem}' for elem in range(1, 14)]
    df_len = pd.DataFrame(columns=cols_len)
    lens_hels = None
    try:
        lens_hels = calc_len_of_hel(pdb_file, ref)
    except KeyError:
        print('KeyError while calculating len of hel')
        pass
    data_lens = [protein_name]
    if lens_hels is not None:
        for lens in lens_hels:
            data_lens.append(lens)
    df_len = df_len.append(pd.Series(data_lens, index=cols_len[0:len(data_lens)]), ignore_index=True)
    return df_len
