"""
Main script for extracting compound/fragment fingerprints.

Inspired by https://github.com/mwinokan/HIPPO
"""
if __name__ == '__main__':
	# Conditional imports only when running as the main script
	from FragFeatures.core.tools import sanitise_mol
	from FragFeatures.utils import timefunction
else:
	from FragFeatures.core.tools import sanitise_mol
	from FragFeatures.utils import timefunction

# TODO: Rename this to ...

import os
import molparse as mp
# from molparse.rdkit.features import FEATURE_FAMILIES, COMPLEMENTARY_FEATURES
from molparse.rdkit import draw_mols

# RDKIT imports
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import MolFromSmiles
import numpy as np
import MDAnalysis as mda
import prolif as plf

import logging
logging.getLogger('MDAnalysis').setLevel(logging.WARNING)
logger = logging.getLogger('FragFeatures') # NOTE: Implement this a bit more


FDEF = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
FEATURE_FAMILIES = FDEF.GetFeatureFamilies()

# COMPLEMENTARY_FEATURES = {
# 	"Donor": "Acceptor",
# 	"Acceptor": "Donor",
# 	"NegIonizable": "PosIonizable",
# 	"PosIonizable": "NegIonizable",
# 	"Aromatic": "Aromatic",
# 	"Aromatic": "PosIonizable",
# 	"PosIonizable": "Aromatic",
# 	"Hydrophobe": "Hydrophobe",
# }
# TODO: Can a frozenset be used here?
COMPLEMENTARY_FEATURES = {
	"Donor": "Acceptor",
	"Acceptor": "Donor",
	"NegIonizable": "PosIonizable",
	"PosIonizable": ["NegIonizable","Aromatic"],
	"Aromatic": "Aromatic",
	"Aromatic": "PosIonizable",
	"Hydrophobe": "Hydrophobe",
}

CUTOFF_PADDING = 0.5

# TODO: Convert the keys to frozensets
FEATURE_PAIR_CUTOFFS = {
	'Donor Acceptor': 3.5 + CUTOFF_PADDING,
	'Acceptor Donor': 3.5 + CUTOFF_PADDING,
	'NegIonizable PosIonizable': 4.5 + CUTOFF_PADDING,
	'PosIonizable NegIonizable': 4.5 + CUTOFF_PADDING,
	'Aromatic PosIonizable': 4.5 + CUTOFF_PADDING,
	'PosIonizable Aromatic': 4.5 + CUTOFF_PADDING,
	'Aromatic Aromatic': 6.0 + CUTOFF_PADDING,
	'Hydrophobe Hydrophobe': 4.5 + CUTOFF_PADDING,
}
# TODO: Evaluate cutoffs and feature selection
# TODO: Allow for easy mofdiciation of feature type selection
# TODO: Consider a an alternative way to define directory structures
class Pose:
	"""
	Represents a pose of a compound in a target.
	"""
	def __init__(self, target_dir, compound_code, interaction_types="hbonds+"):
		self.target_dir = target_dir
		self.compound_code = compound_code
		self.id = compound_code # TODO: for now: what is the id for a pose? Ask Max
		self.protein_path = f'{target_dir}/aligned_files/{compound_code}/{compound_code}_apo.pdb'
		self.mol_path = f'{target_dir}/aligned_files/{compound_code}/{compound_code}_ligand.mol'
		self.smi_path = f'{target_dir}/aligned_files/{compound_code}/{compound_code}_ligand.smi'
		self.protein_system = mp.parse(self.protein_path, verbosity=False).protein_system
		self.feature_families, self.complementary_features, self.ligand_families = self.define_interaction_types(interaction_types)

		# self.complex_system = None
		# self.reference = None
		# self.inspirations = []
		# self.fingerprint = None
		# self.features = self.calculate_features()

		self.sparse_fingerprint = None # TODO: Is this the way to store the fingerprint? Ask Max
		self.sparse_fingerprint_ext = None
		self.fingerprint = None
		self.fingerprint_ext = None

		self._mol = None  # TODO: Why is there an underscore here? Ask Max


	def define_interaction_types(self, selected_interaction_types):
		"""Define the interaction types for the pose by picking feature families"""
		# Check if the selected interaction types are valid
		feature_families = list(FEATURE_FAMILIES)
		complementary_features = dict(COMPLEMENTARY_FEATURES)
		
		valid_interactions = ['hbonds', 'hbonds+', '-aromatic', 'all', ]
		if selected_interaction_types not in valid_interactions:
			raise ValueError(f'Invalid interaction types: {selected_interaction_types}. Valid options are: {valid_interactions}')
		
		# Define the interaction types
		if selected_interaction_types == 'hbonds':
			interaction_types = ['Donor', 'Acceptor']
		elif selected_interaction_types == 'hbonds+':
			interaction_types = ['Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'PosIonizable']
		elif selected_interaction_types == '-aromatic':
			interaction_types = ['Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'PosIonizable', 'Hydrophobe']
		elif selected_interaction_types == 'all':
			interaction_types = list(FEATURE_FAMILIES) # TODO: Test this properly
		else:
			raise NotImplementedError
		
		# Modify the feature families list to include only defined interaction types
		# Allows for duplicates
		updated_feature_families = [f for f in feature_families for _ in range(interaction_types.count(f))]

		# Modify the complementary features dict to include only defined interaction types
		updated_complementary_features = {f: complementary_features[f] for f in updated_feature_families for _ in range(interaction_types.count(f))}

		# Get relevant feature famillies for the ligand
		ligand_families = []
		for sublist in updated_complementary_features.values():
			if isinstance(sublist, str):
				ligand_families.append(sublist)
			else:
				ligand_families.extend(sublist)  # if a list, unpack

		return updated_feature_families, updated_complementary_features, ligand_families

	@timefunction
	def calculate_fingerprint(self, verbose=False):
		"""Calculate the pose's interaction fingerprint"""
		#FIXME: Get this working with the protonoted ligand and protein
		# TODO: Add some error handling, type checking, and logging

		protein_system = self.protein_system
		assert protein_system  # NOTE: What exatly is this doing?

		if not self.mol:
			return

		comp_features = self.features # TODO: What is going on here? Ask Max

		# feature_families = self.feature_families
		ligand_families = self.ligand_families
		compound_features_by_family = {}
		for family in ligand_families: # NOTE: was feature families before
			compound_features_by_family[family] = [f for f in comp_features if f.family == family]

		# Use molparse to get protein features
		protein_features = self.protein_system.get_protein_features()

		sparse_fingerprint = {}
		sparse_fingerprint_ext = {} # extended sparse fingerprint

		chains = protein_system.chain_names
		print(f'\nCompound: {self.compound_code}')
		print(f'Chains: {chains}')

		# 
		for prot_feature in protein_features:

			if prot_feature.res_chain not in chains:
				print(f'Chain {prot_feature.res_chain} not in chains')
				continue

			# Check if the feature is in a user-selected family
			if prot_feature.family not in self.feature_families:
				continue

			# Get the feature family and residue from the protein_system obj
			prot_family = prot_feature.family
			prot_residue = protein_system.get_chain(prot_feature.res_chain).residues[f'n{prot_feature.res_number}']

			if not prot_residue:
				continue

			if prot_residue.name != prot_feature.res_name:
				logger.warning(f'Feature {repr(prot_feature)}')
				continue

			atom_names = [a.name for a in prot_feature.atoms]
			prot_atoms = [prot_residue.get_atom(a) for a in atom_names]

			prot_coords = [a.np_pos for a in prot_atoms if a is not None]
			prot_coord = np.array(np.sum(prot_coords,axis=0)/len(prot_atoms))

			complementary_family = self.complementary_features[prot_family]
			# Check if the family has multiple potential complementary features
			if isinstance(complementary_family, list):
				for i in complementary_family:
					try:
						# print(i)
						complementary_comp_features = compound_features_by_family[i]
						cutoff = FEATURE_PAIR_CUTOFFS[f'{prot_family} {i}']
					except:
						pass
			else:
				complementary_comp_features = compound_features_by_family[complementary_family]
				cutoff = FEATURE_PAIR_CUTOFFS[f'{prot_family} {complementary_family}']

			valid_features = [f for f in complementary_comp_features if np.linalg.norm(f - prot_coord) <= cutoff] # TODO: understand this line better

			# Naming the feature (id)
			if len(prot_feature.atoms) > 1:  # protocol if more than one atom for protein feature
				feature_duck_name = f"{prot_feature.res_chain}_{prot_feature.res_name}_{prot_feature.res_number}_{prot_feature.atoms[0]}" # TODO: implement this properly. 
				if len(valid_features) >= 1:  # print warning if feature is valid
					logger.warning(f"More than one ligand atom associated with protein feature {feature_duck_name}: {prot_feature.atoms}") # TODO: Add compound code here
					# TODO: Put the warning in the metadata
				else:
					pass
			else:
				feature_duck_name = f"{prot_feature.res_chain}_{prot_feature.res_name}_{prot_feature.res_number}_{prot_feature.atoms[0]}"

			# TODO: This needs to be made clearer
			sparse_fingerprint[prot_feature.family_name_number_chain_atoms_str] = (len(valid_features), feature_duck_name) # just the number of ligand atoms + name
			sparse_fingerprint_ext[prot_feature.family_name_number_chain_atoms_str] = (prot_feature, valid_features, feature_duck_name) # extended: includes feature details

		self.sparse_fingerprint = sparse_fingerprint # TODO: Is this the way to store the fingerprint? Ask Max
		self.sparse_fingerprint_ext = sparse_fingerprint_ext

		# More useful fingerprint dict with only values > 0 == interaction present
		dense_fingerprint = {k: v for k, v in sparse_fingerprint.items() if v[0] > 0}
		dense_fingerprint_ext = {k: v for k, v in sparse_fingerprint_ext.items() if len(v[1]) > 0}

		# dict -> {family + protein feature: (ligand atom count, feature_duck_name)}
		self.fingerprint = dense_fingerprint
		# dict -> {family + protein feature: (protein_feature, ligand atoms, feature_duck_name)}
		self.fingerprint_ext = dense_fingerprint_ext

		# List each dictionary key-value pair on a new line for printing
		dense_fingerprint_str = '\n'.join([f'{k}: {v}' for k, v in dense_fingerprint.items()])
		dense_fingerprint_ext_str = '\n'.join([f'{k}: {v}' for k, v in dense_fingerprint_ext.items()])

		print(f'Number of features: {len(dense_fingerprint.items())}\n')

		if verbose:
			# Print the fingerprint
			print(f"\nFingerprint\n-----------\n{dense_fingerprint_str}")
			print(f"\nDense Fingerprint\n-----------------\n{dense_fingerprint_ext_str}\n")


	@timefunction
	def calculate_prolif_fp(self, output_dir, rdkit_protein=False):
		"""
		Calculate the ProLIF fingerprint for a pose.
		"""
		# Add hydrogens to the protein
		from FragFeatures.sanitisation.protein_preparation import ProteinPreparation # FIXME: Move this to the top

		# Add some timing
		import time
		start_time = time.time()


		# FIXME: Uncomment this
		# protein_prep = ProteinPreparation(protein_path=self.protein_path, output_dir=output_dir)
		# protein_prep.prepare_protein()
		# prepared_protein_path = protein_prep.get_prepared_protein_path()
		# # print(f'Prepared protein path: {prepared_protein_path}')

		# TEMPORARY # FIXME: Remove this
		prepared_protein_path = '/Users/nfo24278/Documents/dphil/diamond/DuCK/code/features/prolif_testing/cx0270a_apo_prepared.pdb'

		if rdkit_protein:
			from rdkit import Chem
			rdkit_prot = Chem.MolFromPDBFile(prepared_protein_path, removeHs=False)
			protein_mol = plf.Molecule(rdkit_prot)
		else:
			u = mda.Universe(prepared_protein_path)
			protein_mol = plf.Molecule.from_mda(u)


		# Timing
		mda_time = time.time()
		print(f"\nplf protein: {mda_time - start_time:.4f} seconds")


		print(f"\nNumber of resiudes: {protein_mol.n_residues}")

		# Timing
		prolif_time = time.time()
		print(f"\n`plf.Molecule.from_mda(u)`: {prolif_time - mda_time:.4f} seconds")

		# use default interactions
		# fp = plf.Fingerprint()
		# TODO: Add option to specify interactions
		# fp = plf.Fingerprint(["HBDonor", "HBAcceptor", "CationPi", "Cationic", "Anionic"])
		fp = plf.Fingerprint() # FIXME: just for testing

		# Timing
		prolif_fp__time = time.time()
		print(f"\n`plf.Fingerprint()`: {prolif_fp__time - prolif_time:.4f} seconds")

		ligand_mol = self.ligand_preparation(output_dir=output_dir)

		# Timing
		ligand__time = time.time()
		print(f"\n`self.ligand_preparation(output_dir=output_dir)`: {ligand__time - prolif_fp__time:.4f} seconds")

		ligand = plf.Molecule.from_rdkit(ligand_mol)


		# Timing
		ligand_plf__time = time.time()
		print(f"\n`plf.Molecule.from_rdkit(ligand_mol)`: {ligand_plf__time - ligand__time:.4f} seconds\n")


		# run on your poses
		fp.run_from_iterable([ligand], protein_mol)

		# Timing
		fp_iter_time = time.time()
		print(f"\n`fp.run_from_iterable([ligand], protein_mol)`: {fp_iter_time - ligand_plf__time:.4f} seconds")


		print('\n')
		# fp.ifp[0].ResidueId
		df = fp.to_dataframe()
		print(f"\nPlain dataframe (.T)\n{df.T}") # TODO: Save this to a file

		# NOTE: This might be the best approach
		print(f"\nColumns\n{df.columns.tolist()}")
		# print(dir(df.columns))
		print(f"\nIndex\n{df.index}")


		print('\n\nACCESSING THE DATA...\n')
		interactions = []
		# NOTE: This is the best way to access the data
		for (lig_res, prot_res, interaction) in df.columns:
			prot_res_result = fp.ifp[0][(lig_res, prot_res)]
			if len(prot_res_result) > 1:
				# print(f"\n{prot_res_result}")
				# print(len(prot_res_result))
				# Loop through all the items in the dictionary
				# print('\niterating...')
				for key, value in prot_res_result.items():
					prot_res_result_i = {}
					prot_res_result_i[key] = value
					# print(prot_res_result_i)
					interactions.append(prot_res_result_i)
			else:
				# print(f"\n{prot_res_result}")
				# print(len(prot_res_result))
				interactions.append(prot_res_result)

		print('\n\nInteractions\n')
		# print(interactions)


		self.prolif_fp = interactions  # TODO: Turn this into a property

		test_idx = 1

		test_interaction = self.prolif_fp[test_idx]
		print(test_interaction)
		# Print the value of the dictionary
		for key, value in test_interaction.items():
			print(f"\nKey\n{key}")
			print(value[0])
			props = value[0]
			parent_indices = props['parent_indices']
			print(f"\nParent indices\n{parent_indices}")
			ligand_idx = parent_indices['ligand']
			protein_idx = parent_indices['protein']




		u = mda.Universe(prepared_protein_path)
		prot_atom = u.select_atoms(f'index {protein_idx[0]}')
		print(f"\nProtein atom\n{prot_atom}")
		# print(dir(prot_atom))
		print(
			f"\n"
			f"Protein ChainID: {prot_atom.chainIDs[0]}\n"
			f"Protein ResName: {prot_atom.resnames[0]}\n"
			f"Protein ResID: {prot_atom.resids[0]}\n"
			# f"Protein ResNums: {prot_atom.resnums[0]}\n"
			f"Protein AtomName: {prot_atom.names[0]}\n"
			# f"Protein Residues: {prot_atom.residues[0]}\n"
			)

		# def mol_with_atom_index(mol):
		# 	for atom in mol.GetAtoms():
		# 		print(atom.GetIdx(), atom.GetSymbol(), atom.GetSmarts(), atom.GetNumExplicitHs(), atom.GetFormalCharge())
		# 		# atom.GetNeighbors()
		# 		# print(dir(atom))
		
		# mol_with_atom_index(ligand_mol)


		# # Save the highlighted molecule as image
		highlighted_ligand_mol_filename = os.path.join(output_dir, f'{self.id}_highlighted.png')

		self.draw_highlighted_mol(mol=ligand_mol, atom_indices=list(ligand_idx), filename=highlighted_ligand_mol_filename)


		# NOTE: Maybe something useful in here - like ligand_interactions
		# # Get unique ligand IDs
		# all_ligands = df.columns.get_level_values('ligand').T
		# print(f"\n\n\nLigands\n{all_ligands}")
		# # print(dir(unique_ligands[0]))

		# interaction_types = df.columns.get_level_values('interaction').T
		# print(f"\nInteraction Types\n{interaction_types}")


		# # Get unique protein residues
		# protein_residues = df.columns.get_level_values('protein')
		# protein_residues_list = protein_residues.tolist()
		# print(f"\nProtein residues list\n{protein_residues_list}")

		# # Basically need to combine the ligand and protein residues to get the interactions
		# # Extract resids (parent_indices), distances and angles


		# # Dictionary to store interactions data for each ligand
		# ligand_interactions = {}

		# # Looping through each ligand to extract their interactions
		# for ligand in all_ligands:
		# 	ligand_interactions[ligand] = df.xs(ligand, level='ligand', axis=1).T

		# print(f"\nLigand Interactions\n{ligand_interactions}")

		# Timing
		end_time = time.time()
		print(f"\n`calculate_prolif_fp` executed in {end_time - start_time:.4f} seconds")


	@timefunction
	def ligand_preparation(self, verbose=False, output_dir=None):
		"""
		Prepare the ligand for featrure analysis.
		"""
		m = Chem.MolFromMolFile(self.mol_path)
		# m = AllChem.AssignBondOrdersFromTemplate(m, )
		if verbose:
			print(f"NumAtoms: {m.GetNumAtoms()}")
		# Chem.AllChem.EmbedMolecule(m)

		m = Chem.AddHs(m, addCoords=True)
		if verbose:
			print(f"NumAtomswithHs: {m.GetNumAtoms()}")
		# Embed the 3D structure

		# Write molecule
		if output_dir:
			mol_path = os.path.join(output_dir, f'{self.id}_Hs.mol')
			out = Chem.SDWriter(mol_path)
			out.write(m)
			out.close()

		return m


	@timefunction
	def draw_highlighted_mol(self, mol, atom_indices, filename=None, img_size=(400, 400), save_3d_img=True):
		"""
		Draw the pose's rdkit.Chem.Mol with highlighted atoms.
		"""
		filename3d = filename.replace('.png', '_3d.png')
		bond_indices = []
		for i in range(len(atom_indices)):
			for j in range(i + 1, len(atom_indices)):
				bond = mol.GetBondBetweenAtoms(atom_indices[i], atom_indices[j])
				if bond is not None:
					bond_indices.append(bond.GetIdx())

		ligand_mol2d = Chem.Mol(mol) # Copy the molecule
		ligand_mol2d.Compute2DCoords()

		img = Draw.MolToImage(
			ligand_mol2d,
			highlightAtoms=atom_indices,
			highlightBonds=bond_indices,
			size=img_size
			)

		if filename:
			img.save(filename)
			if save_3d_img:
				img3d = Draw.MolToImage(
					mol,
					highlightAtoms=atom_indices,
					highlightBonds=bond_indices,
					size=img_size
					)
				img3d.save(filename3d)
		else:
			img.show()


	# TODO: Delete?
	def draw_mol3d(self, filename=None):
		"""
		Draw the pose's rdkit.Chem.Mol. This is drawn considering 3D coords.
		"""
		if filename:
			Draw.MolToFile(self.mol, filename)
		else:
			img = Draw.MolToImage(self.mol)
			img.show()


	# TODO: Delete?
	def draw_mol(self, filename=None):
		"""
		Draw the pose's rdkit.Chem.Mol. This is drawn considering 2D coords.
		"""
		if filename:
			m = MolFromSmiles(self.smiles)
			Draw.MolToFile(m, filename)
		else:
			m = MolFromSmiles(self.smiles)
			img = Draw.MolToImage(m, size=(300, 300))
			img.show()


	def draw(self, inspirations=True, protein=False, **kwargs):
		"""Render this pose (and its inspirations)"""
		
		if protein:
			from molparse.py3d import render
			return render(self.complex_system, **kwargs)

		mols = [self.mol]
		if inspirations and self.inspirations:
			mols += [i.mol for i in self.inspirations]

		return draw_mols(mols)


	@property
	def duck_feature_names(self):
		"""Return the pose's protein feature names in a list format suitable for the DUck output"""
		# Combine only the second index of the dict valies into a LIST for DUck
		return [v[1] for v in self.fingerprint.values()]
		# return self.fingerprint_ext


	@property
	def protein_features(self):
		"""Return the pose's protein features as a list"""
		return [v[0] for v in self.fingerprint_ext.values()]


	@property
	def ligand_features(self):
		"""Return the pose's ligand features as a list"""
		return [v[1] for v in self.fingerprint_ext.values()]


	@property
	def smiles(self):
		"""Return the pose's SMILES"""
		with open(self.smi_path, 'r') as f:
			return f.read().strip()


	@property
	def mol(self):  # TODO: compare this to the mol property in pose.py
		"""Returns a pose's rdkit.Chem.Mol"""
		if not self._mol:
			if self.mol_path.endswith('.mol'):
				logger.reading(self.mol_path)
				mol = mp.parse(self.mol_path, verbosity=False)

				if not mol:
					logger.error(f'[{self}] Error computing RDKit Mol from .mol={self.path}')
					raise InvalidMolError

				self.mol = mol

			else:
				raise NotImplementedError

			if not mol:
				logger.error(f'Could not parse {self}.path={self.mol_path}')

		return self._mol


	@mol.setter
	def mol(self, m):
		"""Set the pose's rdkit.Chem.Mol"""
		assert m
		self._mol = sanitise_mol(m)
		# self.db.update(table='pose', id=self.id, key='pose_mol', value=m.ToBinary())


	@property
	def features(self):
		"""Returns the pose's features"""
		return mp.rdkit.features_from_mol(self.mol)
	

	# NOTE: Look into setting a path for when this is all run in a script


	# def calculate_prolif_fp(self):
	# 	"""
	# 	Calculate the PProLIF fingerprint for a pose.
	# 	"""
	# 	# NOTE: Need to add hydrogens first
	# 	# Create an MDAnalysis Universe
	# 	u = mda.Universe(self.protein_path)
	# 	protein_mol = plf.Molecule.from_mda(self.mol)


class InvalidMolError(Exception):
	...


# Test
if __name__ == '__main__':
	from FragFeatures.target_parser import TargetParser
	target = TargetParser('/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac')
	lig = Pose(target.target_dir, 'cx0270a')
	# lig.calculate_fingerprint()
	# print(lig.duck_feature_names)
	# print(target.target_dir)
	# print(lig.protein_path)
	# print(lig.mol_path)
	print('\n')
	lig.calculate_prolif_fp(
		'/Users/nfo24278/Documents/dphil/diamond/DuCK/code/features/Protein_preparation',
		rdkit_protein=True
		)
	
