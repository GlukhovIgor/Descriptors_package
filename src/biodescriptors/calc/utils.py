import yaml
import json
import copy
import logging
from typing import Union

from Bio import PDB

logger = logging.getLogger(__name__)


#_________________FUNCTIONS____________________________#

def load_config(config_file: Union[str, dict]) -> dict:
    """
    Loads file with some type of configuration.
    Config types:
        Descriptor config - contains descriptor settings

    Examples of all configs can be found in /configs.

    Parameters:
    ----------
    config_file: str or dict
        If str is passed, config_file is interpreted as path to a file.
        A file should be of json or yaml format to be correctly loaded.
        If dict is passed, the function returns a deep copy of an object.
    """
    if isinstance(config_file, str):
        if config_file.endswith(".yml") or config_file.endswith(".yaml"):
            try:
                with open(config_file, "r") as stream:
                    return yaml.safe_load(stream)
            except BaseException as e:
                logger.error(f"Failed to read yaml file. Please check format validity at {config_file}.\n{e}")
                raise
        elif config_file.endswith(".json"):
            try:
                with open(config_file, "r") as f:
                    return json.load(f)
            except BaseException as e:
                logger.error(f"Failed to read json file. Please check format validity at {config_file}.\n{e}")
                raise
    elif isinstance(config_file, dict):
        return copy.deepcopy(config_file)
    else:
        logger.error(
            """Invalid format! Please supply a path to .yml, .yaml, or .json file as str.
            Alternatively, provide python dict as an argument."""
        )
    return None


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