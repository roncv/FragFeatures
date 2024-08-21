"""
Main script for parsing through a Fragalysis target's directory. 

Relies on the directory structure of a Fragalysis target.
"""
if __name__ == '__main__':
	# Conditional imports only when running as the main script
	from tools import sanitise_mol
	from FragFeatures.pose import Pose
else:
	from FragFeatures.tools import sanitise_mol

import molparse as mp
from molparse.rdkit.features import FEATURE_FAMILIES, COMPLEMENTARY_FEATURES
from molparse.rdkit import draw_mols
# from molparse.list import NamedList
# from molparse.group import AtomGroup

# from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
import pandas as pd

# from pathlib import Path
import logging
logger = logging.getLogger('FragFeatures')

# TODO: Option to specify location of metadata.csv
class TargetParser:
	"""
	Parse through an XChem target's directory.

	Directory given should be the parent directory containing the target's metadata.csv file.
	"""
	def __init__(self, target_dir):
		self.target_dir = target_dir
		self.metadata = pd.read_csv(f'{target_dir}/metadata.csv') #, index_col='Code') # Code is the unique identifier for each compound

	def get_all_compounds(self):
		"""
		Return all compounds in the target.
		"""
		return self.metadata['Code'].tolist()

	def get_compound(self, code):
		"""
		Return a compound's metadata with the given code.
		"""
		pass


# Test
if __name__ == '__main__':
	target = TargetParser('/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac')
	# print(target.metadata.head())
	lig = Pose(target.target_dir, 'cx0294a')
	# print(lig.mol)
	#lig.draw_mol()
	lig.calculate_fingerprint()
	print(lig.duck_features)
	print(target.target_dir)
	print(lig.protein_path)
	print(lig.mol_path)

	# print(lig.fingerprint)