# FragFeatures

# Welcome to FragFeatures

Tool to extract protein-ligand interaction fingerprints from Fragalysis data.

## Installation

### Conda

It's recommended that you work in a conda environment. To create a new environment and install the required packages, run the following commands:
```{bash}
$ git clone https://github.com/roncv/FragFeatures.git
$ conda create -n fragfeatures -c conda-forge python=3.11
$ conda activate fragfeatures
$ conda install -c conda-forge numpy rdkit pandas pyyaml
$ pip install molparse
$ cd FragFeatures
$ python setup.py install
```
