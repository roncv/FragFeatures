"""
Process the output of a compound's DUck simulation.
"""

import os
import csv

class DUckCompound():
	"""
	Process the output of a compound's DUck simulation.
	"""
	def __init__(
			self,
			experiment_dir,
			compound_id=None,
			output_dir='Analysis',
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


	def validate_simulation(self):
		"""
		Validate that the simulation has run.
		"""
		pass


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
				raise ValueError("No WQB found in the file.")


	def get_wqbs(self):
		"""
		Return a list of wqb values for all simulations of a compound.
		"""
		wqbs = {}
		for feature_dir, feature_name in zip(self.feature_dirs, self.simulation_dirnames):
			wqb = self.get_wqb(feature_dir=feature_dir)
			wqbs[feature_name] = wqb

		return wqbs


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

		return compound_dir


	@property
	def simulation_dirnames(self):
		"""
		Return a list of all feature directories for a compound with simulation data.
		"""
		compound_dir = self.compound_dir
		simulation_dirnames = [d for d in os.listdir(compound_dir) if os.path.isdir(os.path.join(compound_dir, d))]
		# Check for the presence of a simulation directory
		simulation_dirnames = [d for d in simulation_dirnames if os.path.isdir(os.path.join(compound_dir, d, self.duck_sim_dirname))]
		if not simulation_dirnames:
			raise FileNotFoundError(f"No simulation directories found for compound '{self.compound_id}'.")

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
			raise FileNotFoundError(f"No simulation directories found for compound '{self.compound_id}'.")

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
			raise FileNotFoundError(f"No simulation directories found for compound '{self.compound_id}'.")

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
	print(f"\n{duck.get_wqbs()}\n")