# Descriptors_package
This repository contains **Biodescriptors** python package. This tool provides functions for calculating various geometrical descriptors for proteins given their structure as PDB-file and outputting them in a form of Dataframe.

## Structure of a package
Currently, this package contains a single subpackage. Below, the structure of the package is listed.


- **`calc`**\
Subpackage that contains functions regarding structural descriptors calculation. This module includes the following scripts:
  - **`calc_all.py`**\
Contains function to select, calculate and save descriptors in a *.csv* format.
  - **`calc_prot_hel_dist.py`**\
Contains functions for calculation of distances between protein's center of mass and between every helix's center of mass.
  - **`calc_pairwise_sep_dist.py`**\
Contains functions for calculation of separation distances between every helix.
  - **`calc_COM_Calpha_angles.py`**\
Contains functions for calculation of angles between protein's center of mass and alpha carbon atom of every helix.
  - **`calc_len_of_hel.py`**\
Contains functions for calculation of lengths of helices from structure.
  - **`calc_angles_between_hel.py`**\
Contains functions for calculation of angles between all helices in structure.
  - **`calc_COM_clamp.py`**\
Contains functions for calculation of distances between protein's center of mass and every charge clamps.
  - **`calc_charge_clamp_dist.py`**\
Contains functions for calculation of distance between charge clamp residues.
  - **`calc_charge_clamp_angles.py`**\
Contains functions for calculation of angles between charge clamp residues.
  - **`calc_acc_per_hel.py`**\
Contains functions for calculation of solvent-accessibility area per helix, Å2.
  - **`calc_dssp_hel.py`**\
Contains functions for calculation of descriptors that portray how much longer/shorter helices are when compared to a reference.
  - **`calc_sse_content.py`**\
Contains functions for calculation of secondary structure content.
  - **`calc_COM_helix.py`**\
Contains functions for calculation of center of mass for every helix in PDB structure. *Currently under development.*
  - **`calc_COM_for_planes.py`**\
Contains functions for calculation of centers of mass for every "sandwich layer" of VDR structure. *Currently under development.*
  - **`calc_COM_protein.py`**\
Contains functions for calculation of protein's center of mass. 
  - **`calc_plane_angles.py`**\
Contains functions for calculation of angles between every layer. l1, l2, l3 - lists which contain helices numbers for every "sandwich layer" of VDR structure. *Currently under development.*
  - **`constraints.py`**\
Contains some commonly used objects, like table of atomic weights.
  - **`utils.py`**\
Contains several functions which are commonly used within the package.

## Installation

To use this package, clone this repository and enter the following command from the location of this project

  `pip install . -r requirements.txt`

Several descriptors require DSSP module to be installed, this can be done using:

`sudo apt-get install dssp`

## Usage

To see an example of usage of this package, see [`usage_example.ipynb`](usage_example.ipynb)
