"""
Main script for parsing through a Fragalysis target's directory and extract compound/fragment fingerprints.

Inspired by 
"""
if __name__ == '__main__':
	# Conditional imports only when running as the main script
	from parse_target import TargetParser, Pose
	from utils import timeit
else:
	from FragFeatures.parse_target import TargetParser, Pose
	from FragFeatures.utils import timeit # NOTE: Necessary?


# import os
# import sys
import shutil

# import molparse as mp

# import numpy as np
# import pandas as pd

from pathlib import Path
import logging
logger = logging.getLogger('FragFeatures')


class DUckInput():
	"""
	Prepare the input for DUck simulation.
	"""
	def __init__(self,
			  compound_selection: list | str,
			  experiment_name: str,
			  target_dir: str,
			  ):
		self.compound_selection = compound_selection
		self.experiment_name = experiment_name
		self.target_dir = target_dir
		self.target = TargetParser(target_dir)
		self.compound_codes = self.validate_compounds()

	# @timeit
	def validate_compounds(self):
		"""
		Get the validate the given compounds against the target.
		"""
		all_compound_codes = self.target.get_all_compounds()
		if isinstance(self.compound_selection, str):
			if self.compound_selection == 'all':
				self.compound_codes = all_compound_codes
				print(f"Selected compound: {self.compound_codes}")
				return self.compound_codes

			else:
				# Check if the given compound code is in the target
				if self.compound_selection not in all_compound_codes:
					raise ValueError(f"Invalid compound code: {self.compound_selection}")
				else:
					self.compound_codes = [self.compound_selection]
					print(f"Selected compound: {self.compound_codes}")
					return self.compound_codes

		elif isinstance(self.compound_selection, list):
			# Match given list of compound codes to all compound codes
			self.compound_codes = [c for c in self.compound_selection if c in all_compound_codes]
			invalid_compound_codes = [c for c in self.compound_selection if c not in all_compound_codes]
			print(f"Selected compounds: {self.compound_codes}")
			if invalid_compound_codes:
				raise ValueError(f"Invalid compound codes present: {invalid_compound_codes}")

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
	def prepare_experiment(self):
		"""
		Build & prepare a directory for the experiment with all files and inputs.
		"""
		# Create the experiment directory in current working directory
		experiment_dir = Path(self.experiment_name)
		experiment_dir.mkdir(exist_ok=True)

		# Create the compound directories in the experiment directory
		for compound_code in self.compound_codes:
			compound_dir = experiment_dir / compound_code
			compound_dir.mkdir(exist_ok=True)

			# Get the compound's features
			compound = self.get_compound(compound_code)

			# Get the compound's features
			features = compound.duck_features
			# Copy necessary files from a given directory using general copy function
			shutil.copy2(compound.protein_path, compound_dir) # protien pdb
			shutil.copy2(compound.mol_path, compound_dir) # ligand mol

			# Create subdirectories for the features
			for feature in features:
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

	# @timeit
	def generate_duck_input(self, compound_code, feature, protein_pdb_path, ligand_mol_path, output_dir):
		"""
		Generate the .yaml input for a DUck simulation.
		"""
		# TODO: Add more options for the input file
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

	duck_input = DUckInput(compound_selection=['cx0270a', 'cx0281a'], experiment_name='Experiment', target_dir='/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac')
	print(duck_input.compound_codes)
	duck_input.prepare_experiment()

	# Timing
	end_time = time.time()
	print(f"Executed in {end_time - start_time:.4f} seconds")