import numpy as np
import sympy

from biodescriptors.calc.calc_COM_for_planes import _calc_COM_for_planes
from biodescriptors.calc import utils


def _calc_plane_angles(chain, l1, l2, l3):
    """Calculate angles between every layer.
    l1, l2, l3 - lists which contain helices numbers for every "sandwich layer" of VDR structure."""

    COM_l1 = _calc_COM_for_planes(chain, l1)
    COM_l2 = _calc_COM_for_planes(chain, l2)
    COM_l3 = _calc_COM_for_planes(chain, l3)

    # Calculate plane for every layer
    args = []
    for el in COM_l1:
        args.append(sympy.Point3D(el[0]))
    first_layer = sympy.Plane(*args)

    args = []
    for el in COM_l2:
        args.append(sympy.Point3D(el[0]))
    second_layer = sympy.Plane(*args)

    args = []
    for el in COM_l3:
        args.append(sympy.Point3D(el[0]))
    third_layer = sympy.Plane(*args)

    # Calculate and return angle between every layers
    return [np.degrees(sympy.N(first_layer.angle_between(second_layer))),
            np.degrees(sympy.N(first_layer.angle_between(third_layer))),
            np.degrees(sympy.N(second_layer.angle_between(third_layer)))]


def calc_plane_angles(pdb_file, l1, l2, l3):
    """Calculate angles between every layer.
    l1, l2, l3 - lists which contain helices numbers for every "sandwich layer" of VDR structure.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.
    l1: list of ints
        List of helices numbers for sandwich layer 1.
    l2: list of ints
        List of helices numbers for sandwich layer 1.
    l3: list of ints
        List of helices numbers for sandwich layer 1.

    Returns
    -------
    pandas.DataFrame with calculated descriptor.

    """
    _, _, _, chain, _ = utils.get_model_and_structure(pdb_file)
    return _calc_plane_angles(chain, l1, l2, l3)


def plane_angles_to_pandas(pdb_file, l1, l2, l3, protein_name=None):
    """Putting angles between every layer in pandas dataframe.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.
    l1: list of ints
        List of helices numbers for sandwich layer 1.
    l2: list of ints
        List of helices numbers for sandwich layer 1.
    l3: list of ints
        List of helices numbers for sandwich layer 1.
    protein_name: str, default=None
        Protein name to be added to the resulting dataframe.

    Returns
    -------

    """
# Not implemented
    return None
