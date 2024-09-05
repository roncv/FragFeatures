"""
Main script for parsing through a Fragalysis target's directory.

Relies on the directory structure of a Fragalysis target.
"""

import logging

import pandas as pd

logger = logging.getLogger("FragFeatures")


# TODO: Option to specify location of metadata.csv
class TargetParser:
    """
    Parse through an XChem target's directory.

    Directory given should be the parent directory containing the target's metadata.csv file.
    """

    def __init__(self, target_dir, verbose=False, verbose_l2=False):
        self.target_dir = target_dir
        self.verbose = verbose
        self.verbose_l2 = verbose_l2
        if self.verbose:
            print("\n")
            logger.reading(f"fragalysis target directory: {self.target_dir}")
            logger.info("\nParsing target metadata...")

        self.metadata = pd.read_csv(
            f"{target_dir}/metadata.csv"
        )  # , index_col='Code') # Code is the unique identifier for each compound

    def get_all_compounds(self):
        """
        Return all compounds in the target.
        """
        if self.verbose:
            logger.info("\nGetting all compounds in the target...")
        return self.metadata["Code"].tolist()

    def get_compound(self, code):
        """
        Return a compound's metadata with the given code.
        """
        pass


# Test
if __name__ == "__main__":
    from FragFeatures.pose import Pose

    target = TargetParser(
        "/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac"
    )
    print(target.metadata.columns)
    lig = Pose(target.target_dir, "cx0294a")
    lig.calculate_prolif_fp(
        "/Users/nfo24278/Documents/dphil/diamond/DuCK/code/features/prolif_testing"
    )
    print(lig.duck_feature_names_plf)
    print(target.target_dir)
    print(lig.protein_path)
    print(lig.mol_path)
