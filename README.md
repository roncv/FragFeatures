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
$ conda install -c conda-forge numpy rdkit pandas pyyaml
$ pip install molparse
$ cd FragFeatures
$ pip install .
```

## Usage

FragFeatures can be used as a command line tool or as a Python package. The executable provided is `fragfeat`, which is added to your PATH upon installation of the package with pip.

FragFeatures comes with two submodules to extract features from fragments of a target:

```{bash}
$ conda activate fragfeatures

$ fragfeat --help

usage: fragfeat [-h]  ...

Open-source toolkit for extracting fragment features.

options:
  -h, --help    show this help message and exit

Open-source fragment feature extraction toolkit. Choose one of the following actions::
  
    hello       Hello function from FragFeatures (testing).
    prepare-duck
                Prepare an experiment for DUck.
```

The modules accept inputs as command line arguments, yaml functionality is not yet implemented.

### Prepare DUck input

To extract features from a target and prepare a directory for simulation, you can use the `prepare-duck` module. This module requires the following arguments:

```{bash}
$ fragfeat prepare-duck --help

usage: fragfeat prepare-duck [-h] [-c COMPOUND_SELECTION [COMPOUND_SELECTION ...]] [-e EXPERIMENT_NAME] [-t TARGET_DIR]

options:
  -h, --help            show this help message and exit
  -c COMPOUND_SELECTION [COMPOUND_SELECTION ...], --compound-selection COMPOUND_SELECTION [COMPOUND_SELECTION ...]
                        Compound selection for the experiment.
  -e EXPERIMENT_NAME, --experiment-name EXPERIMENT_NAME
                        Name of the experiment.
  -t TARGET_DIR, --target-dir TARGET_DIR
                        Path to the Fragalysis target's directory.
```

Example:
```
$ fragfeat prepare-duck -c <compound_code(s)> -e Experiment -t /path/to/target/dir
```