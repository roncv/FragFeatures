"""
Prepare the input and directories for DUck simulations.
"""
from FragFeatures.target_parser import TargetParser
from FragFeatures.pose import Pose
from FragFeatures.utils import timefunction, dict_to_json # NOTE: Necessary?

import numpy as np
import shutil
import json
from contextlib import redirect_stdout
import io
from pathlib import Path
import pickle

import logging
logger = logging.getLogger('FragFeatures')


class DUckInput():
	"""
	Prepare the input for DUck simulation.
	"""
	def __init__(
			self,
			compound_selection: list | str,
			experiment_name: str,
			target_dir: str,
			verbose: bool = False,
			verbose_l2: bool = False
			):
		self.compound_selection = compound_selection
		self.experiment_name = experiment_name
		self.target_dir = target_dir
		self.verbose = verbose
		self.verbose_l2 = verbose_l2
		self.target = TargetParser(
			target_dir,
			verbose=self.verbose,
			verbose_l2=self.verbose_l2
			)  # Fragalysis
		self.compound_codes = self.validate_compounds()
		# TODO: Add an option to select all compounds - for script execution

	# @timeit
	def validate_compounds(self):
		"""
		Get the validate the given compounds against the target.
		"""
		all_compound_codes = self.target.get_all_compounds()
		if isinstance(self.compound_selection, str):
			if self.compound_selection == 'all':
				self.compound_codes = all_compound_codes
				if self.verbose:
					print(f"\nSelected compounds: \n{self.compound_codes}")
				return self.compound_codes

			else:
				# Check if the given compound code is in the target
				if self.compound_selection not in all_compound_codes:
					raise ValueError(f"Invalid compound code: \n{self.compound_selection}")
				else:
					self.compound_codes = [self.compound_selection]
					if self.verbose:
						print(f"\nSelected compounds: \n{self.compound_codes}")
					return self.compound_codes

		elif isinstance(self.compound_selection, list):
			# FIXME: Temporary fix for selecting all compounds
			if self.compound_selection == ['all']:
				self.compound_codes = all_compound_codes
				if self.verbose:
					print(f"\nSelected compounds: \n{self.compound_codes}")
				return self.compound_codes

			else:
				# Match given list of compound codes to all compound codes
				self.compound_codes = [c for c in self.compound_selection if c in all_compound_codes]
				invalid_compound_codes = [c for c in self.compound_selection if c not in all_compound_codes]
				if self.verbose:
					print(f"\nSelected compounds: \n{self.compound_codes}")
				if invalid_compound_codes:
					if self.verbose:
						print(f"\nInvalid compound codes present: \n{invalid_compound_codes}")
					raise ValueError(f"Invalid compound codes present: \n{invalid_compound_codes}")
				return self.compound_codes

		else:
			raise ValueError(f"Invalid compound selection: {self.compound_selection}")


	# @timeit
	def get_compound(self, compound_code):
		"""
		Return a compound's features and metadata from the given compound code.
		"""
		compound = Pose(self.target_dir, compound_code)
		compound.calculate_fingerprint()

		return compound

	# @timeit
	def get_prolif_compound(self, compound_code, output_dir):
		"""
		Return a compound's features and metadata from the given compound code.
		"""
		compound = Pose(
			self.target_dir,
			compound_code,
			verbose=self.verbose, 
			verbose_l2=self.verbose_l2
			)
		compound.calculate_prolif_fp(
			output_dir=output_dir
			)

		return compound


	def prepare_experiment_prolif(self):
		"""
		Build & prepare a directory for the experiment with all files and inputs.
		"""
		# Create the experiment directory in current working directory
		experiment_dir = Path(self.experiment_name)
		experiment_dir.mkdir(exist_ok=True)

		compound_summaries = {}
		interaction_types_summary = {}
		compound_tally = 0
		feature_tally = 0

		# Create the compound directories in the experiment directory
		for compound_code in self.compound_codes:
			print(f"\nPreparing {compound_code}...\n")
			compound_dir = experiment_dir / compound_code # TODO: Use Path
			compound_dir.mkdir(exist_ok=True)

			# Get the compound's features
			compound = self.get_prolif_compound(
				compound_code,
				output_dir=compound_dir
				)

			# Get the compound's features
			feature_names = compound.duck_feature_names_plf
			feature_objs = compound.plf_feature_objs

			# Copy necessary files from a given directory using general copy function
			# shutil.copy2(compound.protein_path, compound_dir) # protien pdb
			shutil.copy2(compound.mol_path, compound_dir) # ligand mol

			# # Save the compound object to a pickle file using pickle library
			# with open(f'{compound_dir}/{compound_code}_Pose_obj.pkl', 'wb') as f:
			# 	pickle.dump(compound, f)

			compound_summaries[compound_code] = {
				'protein_path': compound.protein_path,
				'prepared_protein_path': compound.prepared_protein_path,
				'ligand_path': compound.mol_path,
				'prepared_ligand_path': compound.prepared_ligand_path,
				'features': feature_names,
				'num_features': len(feature_names)
			}
			compound_tally += 1
			feature_tally += len(feature_names)

			self.generate_compound_mdata_plf(
				compound=compound,
				compound_code=compound_code,
				compound_dir=compound_dir
				)

			compound.get_interactive_ligand_network(
				compound.plf_fp_obj, 
				output_dir=compound_dir)

			# TODO: Include warning of multiple ligand features
			# TODO: Check if the feature directories already exist and contain simulation data
			# FIXME: Figure out a new naming scheme for the feature directories
			for (feature_name, feature_obj) in zip(feature_names, feature_objs):
				# Increase count of interaction types each time a new one is found
				interaction_types_summary[feature_obj.interaction_type] = interaction_types_summary.get(feature_obj.interaction_type, 0) + 1

				feature_dirname = f"{feature_name}_{feature_obj.interaction_type}"
				feature_dir = compound_dir / feature_dirname
				if self.verbose:
					print("...")
					if self.verbose_l2:
						print(f"\nCreating feature directory: {feature_dir}\n")
				feature_dir.mkdir(exist_ok=True)

				# # Save the feature object to a pickle file using pickle library
				# with open(f'{feature_dir}/{feature_name}_PLFeature_obj.pkl', 'wb') as f:
				# 	pickle.dump(compound, f)

				# Generate the input for DUck simulation
				self.generate_duck_input(
					compound_code=compound_code,
					feature=feature_name,
					protein_pdb_path=f'{compound.prepared_protein_path}',
					ligand_mol_path=f'{compound.prepared_ligand_path}',
					output_dir=feature_dir
					)

				self.generate_feat_mdata_plf(
					feature_name=feature_name,
					feature=feature_obj,
					output_dir=feature_dir
					)
				
				# Generate feature-highlighted images of the ligand molecle
				compound.draw_highlighted_mol(
					atom_indices=list(feature_obj.ligand_indices),
					filename=f'{feature_dir}/{compound_code}_{feature_name}_highlighted.png'
					)

		# Generate metadata for the experiment
		# TODO: Create a large tsv with a summary of all the features
		dict_to_json(compound_summaries, f'{experiment_dir}/compound_summaries.json')
		tallies = {
			'num_compounds': compound_tally,
			'num_features': feature_tally
		}
		dict_to_json(tallies, f'{experiment_dir}/tallies.json')
		dict_to_json(interaction_types_summary, f'{experiment_dir}/interaction_types_summary.json')


	# TODO: Add parameters to deal with overwriting existing directories
	# FIXME: Make this morre modular to allow for integration with other tools
	# @timeit
	def prepare_experiment(self):
		"""
		Build & prepare a directory for the experiment with all files and inputs.
		"""
		# Create the experiment directory in current working directory
		experiment_dir = Path(self.experiment_name)
		experiment_dir.mkdir(exist_ok=True)

		compound_summaries = {}
		compound_tally = 0
		feature_tally = 0

		# Create the compound directories in the experiment directory
		for compound_code in self.compound_codes:
			compound_dir = experiment_dir / compound_code
			compound_dir.mkdir(exist_ok=True)

			# Get the compound's features
			compound = self.get_compound(compound_code) # Run calculate_fingerprint()

			# Get the compound's features
			feature_names = compound.duck_feature_names
			expanded_features = compound.fingerprint_ext # Detailed features
			protein_features = compound.protein_features
			ligand_features = compound.ligand_features

			# Copy necessary files from a given directory using general copy function
			shutil.copy2(compound.protein_path, compound_dir) # protien pdb
			shutil.copy2(compound.mol_path, compound_dir) # ligand mol

			# Save the compound object to a pickle file using pickle library
			with open(f'{compound_dir}/{compound_code}_compound.pkl', 'wb') as f:
				pickle.dump(compound, f)

			compound_summaries[compound_code] = {
				'protein_path': compound.protein_path,
				'ligand_path': compound.mol_path,
				'features': feature_names,
				'num_features': len(feature_names)
			}
			compound_tally += 1
			feature_tally += len(feature_names)

			# FIXME: Include ligand smiles in the compound metadata
			self.generate_compound_metadata(compound=compound,
									 compound_code=compound_code,
									 compound_dir=compound_dir)
			# TODO: Include warning of multiple ligand features
			# Create subdirectories for the features
			# TODO: Check if the feature directories already exist and contain simulation data
			# FIXME: Figure out a new naming scheme for the feature directories
			for (feature,feature_e, prot_feat, lig_feat) in zip(feature_names,
													   expanded_features,
													   protein_features,
													   ligand_features):
				feature_dir = compound_dir / feature
				# print(f"Creating feature directory: {feature_dir}")
				feature_dir.mkdir(exist_ok=True)

				# Generate the input for DUck simulation
				self.generate_duck_input(compound_code=compound_code,
							 feature=feature,
							 protein_pdb_path=f'../{compound_code}_apo.pdb',
							 ligand_mol_path=f'../{compound_code}_ligand.mol',
							 output_dir=feature_dir
							 )

				self.generate_feature_metadata(feature=feature,
								   expanded_feature=feature_e,
								   protein_feature=prot_feat,
								   ligand_features=lig_feat,
								   output_dir=feature_dir
								   )

		# Generate metadata for the experiment
		# TODO: Create a large tsv with a summary of all the features
		dict_to_json(compound_summaries, f'{experiment_dir}/compound_summaries.json')
		tallies = {
			'num_compounds': compound_tally,
			'num_features': feature_tally
		}
		dict_to_json(tallies, f'{experiment_dir}/tallies.json')


	def generate_feat_mdata_plf(
			self,
			feature_name,
			feature,
			output_dir
			):
		"""
		Generate metadata for a feature. Protein -> Multiple liand features.
		"""
		# TODO: Ensure these are compatible with pandas
		feature_details = {
			'DUck_feature': feature_name,
			'prot_feature_family': feature.interaction_type,
			'prot_residue_name': feature.prot_atom_resname,
			'prot_res_num': feature.prot_atom_resid,
			'prot_chain': feature.prot_atom_chainID,
			'prot_atom_name': feature.prot_atom_name,
			'prot_atom_names': list(feature.prot_atom_names),
			'prot_atom_idx': feature.prot_idx,
			'prot_atom_idxs': feature.protein_indices,
			'lig_atom_idxs': feature.ligand_indices
			}

		# Write protein_feature_details dict to a json file 
		dict_to_json(feature_details, f'{output_dir}/feature_metadata.json')

		# Ligand features
		ligand_feature_details = {
			'DUck_feature': feature_name,
			'prot_feature_family': feature.interaction_type,
			'atom_nums': feature.ligand_indices,
			'atom_names': feature.ligand_atom_symbols,
			'atom_charges': feature.ligand_atom_charges,
			'interaction_distance': feature.interaction_distance,
			'dha_angle': feature.dha_angle,
			'atom_smarts': feature.ligand_atom_smarts,
			'smarts_substructure': feature.smarts_substructure,
			'explicit_hydrogens': feature.ligand_atom_num_explicit_hs
			}

		dict_to_json(ligand_feature_details, f'{output_dir}/ligand_feature_metadata.json')


	def generate_feature_metadata(self, feature, expanded_feature,
							   protein_feature, ligand_features, output_dir):
		"""
		Generate metadata for a feature. Protein -> Multiple liand features.
		"""
		# Feature metadata
		feature_details = {
			'DUck_feature': feature,
			'protein_feature': expanded_feature
		}

		# Write feature_details dict to a json file in a pythonic way
		with open(f'{output_dir}/feature_metadata.json', 'w') as f:
			json.dump(feature_details, f, indent=4)

		# NOTE: Needs to be ready for multiple protein features
		# TODO: Ensure these are compatible with pandas
		# Protein feature - not complete (process will need to be reversed)
		# Protein feature metadata
		protein_feature_details = {
			'DUck_feature': feature,
			'feature_family': protein_feature.family,
			'residue_name': protein_feature.res_name,
			'res_num': protein_feature.res_number,
			'chain': protein_feature.res_chain,
			'atom_nums': protein_feature.atom_numbers,
			# 'atoms': protein_feature.atoms,
			'position_xyz': [protein_feature.x,
							 protein_feature.y,
							 protein_feature.z]
		}

		# # Write protein _eature_details dict to a json file 
		dict_to_json(protein_feature_details, f'{output_dir}/protein_feature_metadata.json')

		# Ligand features
		for i, ligand_feature in enumerate(ligand_features):
			# Ligand feature metadata
			# Calculate the distances between the protein and ligand atoms
			get_interaction_distance = lambda atom1, atom2: np.linalg.norm(atom1 - atom2) # cartesian distance
			atom2 = protein_feature.position
			distances = [get_interaction_distance(atom1.position, atom2) for atom1 in ligand_feature.atoms]

			# Ligand feature summaries - redirect stdout to a string
			atom_summaries = []
			for atom in ligand_feature.atoms:
				with redirect_stdout(io.StringIO()) as f:
					atom.summary()
				s = f.getvalue()
				atom_summaries.append(s)

			ligand_feature_details = {
				'DUck_feature': feature,
				'ligand_feature': ligand_feature.__str__(),
				'feature_family': ligand_feature.family,
				'atom_nums': ligand_feature.atom_numbers,
				'atoms_names': [atom.name for atom in ligand_feature.atoms],
				'atoms_xyz': [[atom.x, atom.y, atom.z] for atom in ligand_feature.atoms],
				'atom_pdb_idxs': [atom.pdb_index for atom in ligand_feature.atoms],
				'interaction_distances': distances,
				'atom_summaries': atom_summaries,
				'feature_xyz': [ligand_feature.x,
								 ligand_feature.y,
								 ligand_feature.z]
			}

			# Write ligand_feature_details dict to a json file in a pythonic way
			dict_to_json(ligand_feature_details, f'{output_dir}/ligand_feature_metadata_{i}.json')


	def generate_compound_mdata_plf(self, compound, compound_code, compound_dir):
		"""
		Generate metadata for the compounds.
		"""
		compound_metadata = {
			'compound_code': compound_code,
			'protein': compound.protein_path,
			'prepared_protein': compound.prepared_protein_path,
			'ligand': compound.mol_path,
			'prepared_ligand': compound.prepared_ligand_path,
			'SMILES': compound.smiles,
			'features': compound.duck_feature_names_plf,
			'num_features': compound.number_of_interactions_plf
		}

		# Write compound_metadata dict to a json file
		dict_to_json(compound_metadata, f'{compound_dir}/{compound_code}_metadata.json')

		# Save an RDKit image of the ligand
		compound.draw_mol3d(f'{compound_dir}/{compound_code}_3d.png')
		compound.draw_mol(f'{compound_dir}/{compound_code}_2d.png')


	# FIXME: This needs to be reverted back to the original version
	def generate_compound_metadata(self, compound, compound_code, compound_dir):
		"""
		Generate metadata for the compounds.
		"""
		compound_metadata = {
			'compound_code': compound_code,
			'protein': compound.protein_path,
			'prepared_protein': compound.prepared_protein_path,
			'ligand': compound.mol_path,
			'prepared_ligand': compound.prepared_ligand_path,
			'SMILES': compound.smiles,
			'features': compound.duck_feature_names_plf,
			'num_features': compound.number_of_interactions_plf
		}

		# Write compound_metadata dict to a json file
		dict_to_json(compound_metadata, f'{compound_dir}/{compound_code}_metadata.json')

		# Save an RDKit image of the ligand
		compound.draw_mol3d(f'{compound_dir}/{compound_code}_3d.png')
		compound.draw_mol(f'{compound_dir}/{compound_code}_2d.png')


	# @timeit
	def generate_duck_input(self, compound_code, feature, protein_pdb_path,
						 ligand_mol_path, output_dir):
		"""
		Generate the .yaml input for a DUck simulation.
		"""
		# TODO: Add more options for the input file
		# TODO: Work on flexibillity and linking to the openduck repo
		# TEST: Test with openduck repo - v. quick simulation - or even just setup / prep?
		input_file = f"""# DuCK input file for {feature} feature from {compound_code}\n
# Main Arguments
interaction : {feature}
receptor_pdb : {protein_pdb_path}
ligand_mol : {ligand_mol_path}
gpu_id : 0

# Chunking Arguments
do_chunk : True
cutoff : 10
ignore_buffers : False

# Preparation Arguments
small_molecule_forcefield : smirnoff
protein_forcefield : amber14-all
water_model : tip3p
ionic_strength : 0.05
waters_to_retain : waters_to_retain.pdb
solvent_buffer_distance : 15
force_constant_eq : 1

# Production Arguments
smd_cycles : 20
md_length : 0.5
wqb_threshold : 6
init_velocities : 0.00001
init_distance : 2.5
fix_ligand : True"""

		# Write the input file
		with open(f'{output_dir}/input.yaml', 'w') as f:
			f.write(input_file)



if __name__ == '__main__':
	# Testing
	# Add some timing
	import time
	start_time = time.time()

	duck_input = DUckInput(
		# compound_selection=['cx0270a', 'cx0281a', 'cx0756c', 'cx0858a'],
		compound_selection=['cx1091b'], # issues with termini
		# compound_selection=['cx0270a'],
		# compound_selection='all',
		experiment_name='Experiment2',
		target_dir='/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac',
		verbose=True,
		verbose_l2=True
		)
	# print(duck_input.compound_codes)
	# duck_input.prepare_experiment()
	duck_input.prepare_experiment_prolif()

	# Timing
	end_time = time.time()
	print(f"Executed in {end_time - start_time:.4f} seconds")