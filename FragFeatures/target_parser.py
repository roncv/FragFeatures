"""
Main script for parsing through a Fragalysis target's directory and extract compound/fragment fingerprints.

Inspired by 
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


# CUTOFF_PADDING = 0.5

# FEATURE_PAIR_CUTOFFS = {
# 	'Donor Acceptor': 3.5 + CUTOFF_PADDING,
# 	'Acceptor Donor': 3.5 + CUTOFF_PADDING,
# 	'NegIonizable PosIonizable': 4.5 + CUTOFF_PADDING,
# 	'PosIonizable NegIonizable': 4.5 + CUTOFF_PADDING,
# 	'Aromatic PosIonizable': 4.5 + CUTOFF_PADDING,
# 	'PosIonizable Aromatic': 4.5 + CUTOFF_PADDING,
# 	'Aromatic Aromatic': 6.0 + CUTOFF_PADDING,
# 	'Hydrophobe Hydrophobe': 4.5 + CUTOFF_PADDING,
# }
# # TODO: Evaluate cutoffs and feature selection
# # TODO: Allow for easy mofdiciation of feature type selection

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

# # TODO: Split Pose and TargetParser into separate files
# class Pose:
# 	"""
# 	Represents a pose of a compound in a target.
# 	"""
# 	def __init__(self, target_dir, compound_code):
# 		self.target_dir = target_dir
# 		self.compound_code = compound_code
# 		self.id = compound_code # TODO: for now: what is the id for a pose? Ask Max
# 		self.protein_path = f'{target_dir}/aligned_files/{compound_code}/{compound_code}_apo.pdb'
# 		self.mol_path = f'{target_dir}/aligned_files/{compound_code}/{compound_code}_ligand.mol'
# 		#self.mol = sanitise_mol(Chem.MolFromMolFile(self.mol_path))
# 		# self.mol = Chem.MolFromMolFile(self.mol_path)
# 		self.protein_system = mp.parse(self.protein_path, verbosity=False).protein_system
# 		# self.complex_system = None
# 		# self.reference = None
# 		# self.inspirations = []
# 		# self.fingerprint = None
# 		# self.features = self.calculate_features()

# 		self.sparse_fingerprint = None # TODO: Is this the way to store the fingerprint? Ask Max
# 		self.sparse_fingerprint_ext = None
# 		self.fingerprint = None
# 		self.fingerprint_ext = None

# 		self._mol = None  # TODO: Why is there an underscore here? Ask Max


# 	def calculate_fingerprint(self, verbose=False):
# 		"""Calculate the pose's interaction fingerprint"""
# 		# if self.path.endswith('.pdb'):

# 		protein_system = self.protein_system
# 		# if not self.protein_system:
# 		# 	# logger.reading(self.path)
# 		# 	protein_system = mp.parse(self.path, verbosity=False).protein_system

# 		# elif self.path.endswith('.mol') and self.reference:

# 		# 	logger.debug('fingerprint from .mol and reference pose')
# 		# 	protein_system = self.reference.protein_system

# 		# else:

# 		# 	logger.debug('calculate_fingerprint()')
# 		# 	raise NotImplementedError

# 		assert protein_system

# 		if not self.mol:
# 			return
		
# 		# print(f'\nMol: {self.mol}')

# 		comp_features = self.features # TODO: What is going on here? Ask Max

# 		comp_features_by_family = {}
# 		for family in FEATURE_FAMILIES:
# 			comp_features_by_family[family] = [f for f in comp_features if f.family == family]

# 		# protein_features = self.target.features
# 		# if not protein_features:

# 		# Use molparse to get protein features
# 		protein_features = self.protein_system.get_protein_features()
# 		# print(protein_features)

# 		sparse_fingerprint = {}
# 		sparse_fingerprint_ext = {} # extended sparse fingerprint

# 		chains = protein_system.chain_names
# 		# residues = protein_system.res_names
# 		# residues_2 = protein_system.residue_names # does the same thing using .append instead of list comprehension 
# 		print(f'\nCompound: {self.compound_code}')
# 		print(f'Chains: {chains}')
# 		# print(f'Residues: {residues}')
# 		# print(f'Residues 2: {residues_2}')

# 		print(f'Number of features: {len(protein_features)}\n')

# 		for prot_feature in protein_features:
# 			# print(dir(f))

# 			if prot_feature.res_chain not in chains:
# 				print(f'Chain {prot_feature.res_chain} not in chains')
# 				continue

# 			# Get the feature family and residue from the protein_system obj
# 			prot_family = prot_feature.family
# 			prot_residue = protein_system.get_chain(prot_feature.res_chain).residues[f'n{prot_feature.res_number}']
# 			# print(f'{prot_residue} {prot_family}')

# 			if not prot_residue:
# 				continue

# 			# # if prot_feature.residue_number == 77:
# 			# # 	logger.debug(repr(prot_feature))

# 			if prot_residue.name != prot_feature.res_name:
# 				logger.warning(f'Feature {repr(prot_feature)}')
# 				continue

# 			atom_names = [a.name for a in prot_feature.atoms]
# 			prot_atoms = [prot_residue.get_atom(a) for a in atom_names]
# 			# print(prot_atoms)
# 			# print(atom_names)

# 			prot_coords = [a.np_pos for a in prot_atoms if a is not None]
# 			prot_coord = np.array(np.sum(prot_coords,axis=0)/len(prot_atoms))

# 			complementary_family = COMPLEMENTARY_FEATURES[prot_family]
# 			complementary_comp_features = comp_features_by_family[complementary_family]
# 			cutoff = FEATURE_PAIR_CUTOFFS[f'{prot_family} {complementary_family}']

# 			valid_features = [f for f in complementary_comp_features if np.linalg.norm(f - prot_coord) <= cutoff] # TODO: understand this line better

# 			# Naming the feature (id)
# 			if len(prot_feature.atoms) > 1:  # protocol if more than one atom for protein feature
# 				feature_duck_name = f"{prot_feature.res_chain}_{prot_feature.res_name}_{prot_feature.res_number}_{prot_feature.atoms[0]}" # TODO: implement this properly. 
# 				if len(valid_features) >= 1:  # print warning if feature is valid
# 					logger.warning(f"More than one atom for pritein feature for {feature_duck_name}: {prot_feature.atoms}") # TODO: Add compound code here
# 				else:
# 					pass
# 			else:
# 				feature_duck_name = f"{prot_feature.res_chain}_{prot_feature.res_name}_{prot_feature.res_number}_{prot_feature.atoms[0]}"

# 			sparse_fingerprint[prot_feature.family_name_number_chain_atoms_str] = (len(valid_features), feature_duck_name) # just the number of ligand features
# 			sparse_fingerprint_ext[prot_feature.family_name_number_chain_atoms_str] = (valid_features, feature_duck_name) # extended: includes feature details
# 			# print(fingerprint)
# 			# print(prot_feature.family_name_number_chain_atoms_str)

# 		self.sparse_fingerprint = sparse_fingerprint # TODO: Is this the way to store the fingerprint? Ask Max
# 		self.sparse_fingerprint_ext = sparse_fingerprint_ext

# 		# More useful fingerprint dict with only values > 0
# 		dense_fingerprint = {k: v for k, v in sparse_fingerprint.items() if v[0] > 0}
# 		dense_fingerprint_ext = {k: v for k, v in sparse_fingerprint_ext.items() if len(v[0]) > 0}

# 		self.fingerprint = dense_fingerprint
# 		self.fingerprint_ext = dense_fingerprint_ext

# 		# List each dictionary key-value pair on a new line for printing
# 		dense_fingerprint_str = '\n'.join([f'{k}: {v}' for k, v in dense_fingerprint.items()])
# 		dense_fingerprint_ext_str = '\n'.join([f'{k}: {v}' for k, v in dense_fingerprint_ext.items()])

# 		if verbose:
# 			# Print the fingerprint
# 			print(f"\nFingerprint\n-----------\n{dense_fingerprint_str}")
# 			print(f"\nDense Fingerprint\n-----------------\n{dense_fingerprint_ext_str}\n")


# 	def draw(self, inspirations=True, protein=False, **kwargs):
# 		"""Render this pose (and its inspirations)"""
		
# 		if protein:
# 			from molparse.py3d import render
# 			return render(self.complex_system, **kwargs)

# 		mols = [self.mol]
# 		if inspirations and self.inspirations:
# 			mols += [i.mol for i in self.inspirations]

# 		return draw_mols(mols)


# 	@property
# 	def duck_features(self):
# 		"""Return the pose's features in a format suitable for the DUck database"""
# 		# return self.fingerprint
# 		# Combine only the second index of the dict valies into a LIST for DUck
# 		return [v[1] for v in self.fingerprint.values()]
# 		# return self.fingerprint_ext


# 	@property
# 	def mol(self):  # TODO: compare this to the mol property in pose.py
# 		"""Returns a pose's rdkit.Chem.Mol"""
# 		if not self._mol:
# 			if self.mol_path.endswith('.mol'):
# 				logger.reading(self.mol_path)
# 				mol = mp.parse(self.mol_path, verbosity=False)

# 				if not mol:
# 					logger.error(f'[{self}] Error computing RDKit Mol from .mol={self.path}')
# 					raise InvalidMolError

# 				self.mol = mol

# 			else:
# 				raise NotImplementedError

# 			if not mol:
# 				logger.error(f'Could not parse {self}.path={self.mol_path}')

# 		return self._mol


# 	@mol.setter
# 	def mol(self, m):
# 		"""Set the pose's rdkit.Chem.Mol"""
# 		assert m
# 		self._mol = sanitise_mol(m)
# 		# self.db.update(table='pose', id=self.id, key='pose_mol', value=m.ToBinary())


# 	@property
# 	def features(self):
# 		"""Returns the pose's features"""
# 		return mp.rdkit.features_from_mol(self.mol)
	
# 	def draw_mol(self):
# 		"""Draw the pose's rdkit.Chem.Mol"""
		
# 		img = Draw.MolToImage(self.mol)
# 		img.show()

# 		# TODO: Add option to save this image to a file
# 	# NOTE: Look into setting a path for when this is all run in a script


# class InvalidMolError(Exception):
# 	...

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