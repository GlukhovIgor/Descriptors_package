'''
This repository contains Biodescriptors python package. 
This tool provides the set of functions for calculating various geometrical descriptors for proteins given their structure as a PDB-file and outputting them in a form of a dataframe or CSV-table.

The structure of the package and detailed description of each function can be found in the current documentation (all the modules listed in the left part of this page).
Currently our package contains only one module named `calc`.

To use this package, clone this repository and enter the following command from the location of this project

`pip install . -r requirements.txt`

Several descriptors require DSSP module to be installed, this can be done using:

`sudo apt-get install dssp`

Additionally one may need KPAX software for providing user settings for protein helices borders: http://kpax.loria.fr/

For developers -

`pip install -e . -r requirements_dev.txt`

To run the calculations one needs to have the data in PDB format and its secondary structure annotation in the list format. 
Charge clamp residues are optional, one can set up a personal configuration file to turn off the calculation of these descriptors.

Detailed examples of usage and the following analysis are located at `usage_example.ipynb` and `analysis_example.ipynb`.

'''

from biodescriptors import calc

__all__ = ["calc"]
