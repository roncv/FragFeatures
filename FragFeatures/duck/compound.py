"""
Process the output of a compound's DUck simulation.
"""
if __name__ == '__main__':
	# Conditional imports only when running as the main script
	from FragFeatures.utils import timeit, dict_to_json
	from FragFeatures.duck.feature import DUckFeature
else:
	from FragFeatures.utils import timeit, dict_to_json
	from FragFeatures.duck.feature import DUckFeature

import os
import csv
import pandas as pd

import logging
logger = logging.getLogger('FragFeatures') # NOTE: Implement this a bit more

class DUckCompound():
	"""
	Process the output of a compound's DUck simulation.
	"""
	def __init__(
			self,
			experiment_dir,
			compound_id=None,
			output_dir='analysis',
			wqb_filename='wqb.txt',
			duck_sim_dirname='duck_runs'
			):
		"""
		Initialise the DUckCompound object.

		compound_id should be the directory name of the compound in the experiment directory.
		"""
		self.experiment_dir = experiment_dir
		self.compound_id = compound_id
		self.wqb_filename = wqb_filename
		self.duck_sim_dirname = duck_sim_dirname
		self.output_dir = output_dir
		# self.compound_dir = None

		# Create the analysis directory for the compound
		self.create_analysis_dir()


	def validate_simulation(self):
		"""
		Validate that the simulation has run.
		"""
		pass


	def create_analysis_dir(self):
		"""
		Create the analysis directory for the compound if it doesn't exist.
		"""
		analysis_dir = os.path.join(self.compound_dir, self.output_dir)
		if not os.path.isdir(analysis_dir):
			os.mkdir(analysis_dir)
		
		self.analysis_dir = analysis_dir


	def get_wqb(self, feature_dir):
		"""
		Return the wqb value of a DUck simulation.
		"""
		# Get the wqb value from the wqb.txt file
		wqb_file = os.path.join(feature_dir, self.wqb_filename)
		with open(wqb_file, newline='') as tsvfile:
			reader = csv.DictReader(tsvfile, delimiter='\t')
			wqb_column_header = reader.fieldnames[1]
			if wqb_column_header != 'WQB':
				print(f"Warning: Second column header is '{wqb_column_header}' instead of 'WQB'.")
			row = next(reader, None)  # read the first (and only) data row
			if row is not None:
				wqb = float(row[wqb_column_header])

				return wqb
			else:
				logger.warning("No WQB found in the file.")
				# raise ValueError("No WQB found in the file.")


	def get_summaries(self):
		"""
		Return a list of wqb values for all simulations of a compound.
		"""
		summaries = []
		# Check that the compound has simulations
		if not self.simulation_dirnames:
			logger.warning(f"No simulations found for compound '{self.compound_id}'.")
			return summaries

		for feature_dir, feature_name in zip(self.feature_dirs, self.simulation_dirnames):
			wqb = self.get_wqb(feature_dir=feature_dir)
			duck_feature = DUckFeature(feature_dir)
			duck_feature.get_feature_metadata()
			# print(duck_feature.feature_dir)
			summaries.append({
				'compound_id': self.compound_id,
				'Feature Name': feature_name,
				'WQB': wqb,
				'protein_feature_family': duck_feature.protein_feature_family,
				'ligand_feature_families': duck_feature.ligand_feature_families,
				'protein_residue': f"{duck_feature.residue_name}{duck_feature.res_num}",
				'protein_chain': duck_feature.res_chain,
				'ligand_atom_names': duck_feature.ligand_atom_names,
				'interaction_types': duck_feature.interaction_types,
				'mixed': duck_feature.mixed,
				'interaction_type': duck_feature.interaction_type
			})

		# Save the summaries to a json file in the analysis directory
		dict_to_json(summaries, os.path.join(self.analysis_dir, 'wqb_summary.json'))

		return summaries


	def plot_simulation_overview(self):
		"""
		Plot the simulation overview. 
		"""
		pass


	def get_simulation_overview(self):
		"""
		Get the summary of a similation. Includes plots, wqb values, etc.
		"""
		pass


	@property
	def compound_dir(self):
		"""
		Get and validate a compounds ID/directory.
		"""
		compound_dir = os.path.join(self.experiment_dir, self.compound_id)
		# Check that the compound_id directory exists
		if not os.path.isdir(compound_dir):
			raise FileNotFoundError(f"Compound directory '{compound_dir}' does not exist.")
			# logger.warning(f"Compound directory '{compound_dir}' does not exist.")

		return compound_dir


	@property
	def simulation_dirnames(self):
		"""
		Return a list of all feature directories for a compound with simulation data.
		"""
		# TODO: Differentiate between feature and simulation directories
		compound_dir = self.compound_dir
		simulation_dirnames = [d for d in os.listdir(compound_dir) if os.path.isdir(os.path.join(compound_dir, d))]
		# Check for the presence of a simulation directory
		simulation_dirnames = [d for d in simulation_dirnames if os.path.isdir(os.path.join(compound_dir, d, self.duck_sim_dirname))]
		if not simulation_dirnames:
			logger.warning(f"No simulation directories found for compound '{self.compound_id}'.")

		return simulation_dirnames


	@property
	def simulation_dirs(self):
		"""
		Return a list of all simulation directories for a compound.
		"""
		compound_dir = self.compound_dir
		simulation_dirnames = [d for d in os.listdir(compound_dir) if os.path.isdir(os.path.join(compound_dir, d))]
		# Check for the presence of a simulation directory
		simulation_dirnames = [d for d in simulation_dirnames if os.path.isdir(os.path.join(compound_dir, d, self.duck_sim_dirname))]
		simulation_dirs = [os.path.join(compound_dir, d, self.duck_sim_dirname) for d in simulation_dirnames]  # full paths
		if not simulation_dirs:
			logger.warning(f"No simulation directories found for compound '{self.compound_id}'.")

		return simulation_dirs


	@property
	def feature_dirs(self):
		"""
		Return a list of all simulation directories for a compound.
		"""
		compound_dir = self.compound_dir
		simulation_dirnames = [d for d in os.listdir(compound_dir) if os.path.isdir(os.path.join(compound_dir, d))]
		# Check for the presence of a simulation directory
		simulation_dirnames = [d for d in simulation_dirnames if os.path.isdir(os.path.join(compound_dir, d, self.duck_sim_dirname))]
		feature_dirs = [os.path.join(compound_dir, d) for d in simulation_dirnames]  # full paths
		if not feature_dirs:
			logger.warning(f"No simulation directories found for compound '{self.compound_id}'.")

		return feature_dirs




# Test
if __name__ == '__main__':
	duck = DUckCompound(
		experiment_dir='/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac_simulations/Experiment',
		compound_id='cx0270a',
		wqb_filename='wqb.txt'
		)
	# print(duck.compound_dir)
	# print(duck.simulation_dirs)
	# print(duck.simulation_dirnames)
	# print(duck.feature_dirs)
	# print('')
	print(f"\n{duck.get_summaries()}\n")

	import pandas as pd
	df = pd.read_json('/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac_simulations/Experiment/cx0270a/analysis/wqb_summary.json')

	# Display the DataFrame
	print(df)