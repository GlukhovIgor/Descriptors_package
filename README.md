# Descriptors_package
This repository contains **Biodescriptors** python package. This tool provides functions for calculating various geometrical descriptors for proteins given their structure as a PDB-file and outputting them in a form of a dataframe or CSV-table.


## Installation

To use this package, clone this repository and enter the following command from the location of this project

  `pip install . -r requirements.txt`

Several descriptors require DSSP module to be installed, this can be done using:

`sudo apt-get install dssp`

Additionally one may need KPAX software for providing user settings for protein helices borders: http://kpax.loria.fr/ 

The structure of the package and detailed description of each function could be found here: https://glukhovigor.github.io/Descriptors_package/biodescriptors.html


## Usage

To see an example of usage of this package, see [`usage_example.ipynb`](usage_example.ipynb)

To see an example of analysis that can be performed using the output of this package, see [`analysis_example.ipynb`](usage_example.ipynb)


## Authors
- Pats Karina, PhD student, ITMO University
- Glukhov Igor, MSc student, ITMO University
- Vinogradova Elisaveta, MSc student, ITMO University


## Acknowledgements
- Dr. Molnar Ferdinand, Associate Professor, Nazarbayev University
- Dr. Marie-Dominique Devinges, Lorraine Research Laboratory in Computer Science and its Applications, University of Lorraine
- Petrosyan Stepan, Bioinformatics Institute student (2019/2020)
- Mamayeva Maria, Bioinformatics Institute student (2019/2020)
