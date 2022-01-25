from Bio import PDB


def get_model_and_structure(pdb_file):
    """
    Initialize PDB structure.
    
    Parameters:
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.

    Returns:
    -------
    tuple: parser, structure, model, chain, atom_structure.

    """
    p = PDB.PDBParser(QUIET=True)
    structure = p.get_structure('protein', pdb_file)
    model = structure[0]
    chain = model['A']
    atom_struct = structure.get_atoms()
    return p, structure, model, chain, atom_struct


def getResidues(dssp):
    """TODO: docstring"""
    residues = []
    numbers = []
    for num, val in enumerate(dssp.keys()):
        numbers.append(num)
        residues.append(val[1][1])
    return residues, numbers


def getNum(n, res_num):
    """TODO: docstring"""
    res, num = res_num
    idx = res.index(n)
    return num[idx]


def getRes(n, res_num):
    """TODO: docstring"""
    res, num = res_num
    idx = num.index(n)
    return res[idx]