# Descriptors_package
This repository contains **Biodescriptors** python package. This tool provides functions for calculating various geometrical descriptors for proteins given their structure and outputting them in a form of Dataframe.

## Structure of a package
- **`calc`**
  - **`calc_all.py`**
Contains function to select, calculate and save descriptors in a *.csv* format.
  - **`calc_{descriptor}.py`**
Multiple files. Each individual file contains calculation logic for its descriptor and a helper function to save its output to a dataframe. 
  - **`constraints.py`**
Contains some commonly used objects, like table of atomic weights.
  - **`utils.py`**
Contains several functions which are commonly used within the package.

## Installation

To use this package, clone this repository and enter the following command from the location of this project

  `pip install . -r requirements.txt`

Several descriptors require DLSS module to be installed, this can be done using:

`sudo apt-get install dlss`
