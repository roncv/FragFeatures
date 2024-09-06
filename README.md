# FragFeatures

# Welcome to FragFeatures

Tool to extract protein-ligand interaction fingerprints from a Fragalysis target.

## Installation

### Conda

It's recommended that you work in a conda environment. To create a new environment and install the required packages, run the following commands:
```{bash}
$ git clone https://github.com/roncv/FragFeatures.git
$ conda create -n fragfeatures -c conda-forge python=3.11
$ conda activate fragfeatures
$ conda install -c conda-forge numpy rdkit pandas openmm
$ pip install molparse
$ cd FragFeatures
$ pip install .
```

## Usage

FragFeatures can be used as a command line tool or as a Python package. The executable provided is `fragfeat`, which is added to your PATH upon installation of the package with pip.

FragFeatures comes with two submodules:

```{bash}
$ conda activate fragfeatures

$ fragfeat --help

usage: fragfeat [-h]                                    ...

Open-source toolkit for extracting fragment features from protein-
ligand complexes.

options:
 -h, --help         show this help message and exit

Subcommands:
                                   
  prepare-duck      Prepare an experiment for DUck.
  summarise-duck    Summarise the output of a DUck experiment.
```

The modules accept inputs as command line arguments, yaml functionality is yet to be implemented.

### Prepare DUck input

To extract features from a target and prepare a directory for simulation, you can use the `prepare-duck` module. This module requires the following arguments:

```{bash}
$ fragfeat prepare-duck --help

usage: fragfeat prepare-duck [-h]
                             [-c COMPOUND_SELECTION [COMPOUND_SELECTION ...]]
                             [-e EXPERIMENT_NAME] [-t TARGET_DIR]
                             [-a] [-v]

options:
  -h, --help            show this help message and exit
  -c COMPOUND_SELECTION [COMPOUND_SELECTION ...], --compound-selection COMPOUND_SELECTION [COMPOUND_SELECTION ...]
                        Compound selection for the experiment.
  -e EXPERIMENT_NAME, --experiment-name EXPERIMENT_NAME
                        Name of the experiment.
  -t TARGET_DIR, --target-dir TARGET_DIR
                        Path to the Fragalysis target's directory.
  -a, --all-compounds   Select all available compounds from a
                        target.
  -v, --verbose         Increase verbosity. Use -v for verbose and
                        -vv for additional verbosity.
```

Example:
```
$ fragfeat prepare-duck -c <compound_code(s)> -e Experiment -t /path/to/target/dir -v
```