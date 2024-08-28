"""
Process the output directory containing results of DUck simulations for multiple compounds.
"""
from FragFeatures.duck.feature import DUckFeature
from FragFeatures.duck.compound import DUckCompound
from FragFeatures.utils import timeit, dict_to_json

import os


class DUckExperiment():
	"""
	Process the output directory containing results of DUck simulations.
	"""
	def __init__(
			self,
			experiment_dir,
			output_dir='analysis',
			wqb_filename='wqb.txt'):
		self.experiment_dir = experiment_dir
		self.wqb_filename = wqb_filename
		self.output_dir = output_dir

		# Create the analysis directory for the experiment
		self.create_analysis_dir()


	def get_simulations(self):
		"""
		Return a list of all simulations in the experiment.
		"""
		pass


	def create_analysis_dir(self):
		"""
		Create the analysis directory for the compound if it doesn't exist.
		"""
		analysis_dir = os.path.join(self.experiment_dir, self.output_dir)
		if not os.path.isdir(analysis_dir):
			os.mkdir(analysis_dir)
		
		self.analysis_dir = analysis_dir


	def summarise_experiment(self):
		"""
		Return a json/csv summary of the experiment.
		"""
		# Combine all the compound summaries into an experiment summary
		experiment_summary = []
		for compound_dir in self.compound_dirs:
			compound_summary = self.get_compound_summary(compound_dir)
			experiment_summary.extend(compound_summary)

		# Write the experiment summary to a json file
		output_file = os.path.join(self.analysis_dir, 'experiment_summary.json')
		dict_to_json(experiment_summary, output_file)
		return experiment_summary


	def get_compound_summary(self, compound_dir):
		"""
		Return a json/csv summary for a specific compound.
		"""
		compound = DUckCompound(
			experiment_dir=self.experiment_dir,
			compound_id=compound_dir,
			wqb_filename=self.wqb_filename
		)
		return compound.get_summaries()


	@property
	def compound_dirs(self):
		"""
		Return a list of all "compound" dirs in the experiment dir.
		"""
		# Get all compound directories in the experiment directory
		# TODO: Tighten this up
		compound_dirs = [dir for dir in os.listdir(self.experiment_dir) if os.path.isdir(os.path.join(self.experiment_dir, dir))]
		if not compound_dirs:
			raise FileNotFoundError("No directories found in the experiment directory.")
		# TODO: Glob for the .pdb files
		# Check that the compound directories contain a .pdb file
		compound_dirs = [dir for dir in compound_dirs if os.path.isfile(os.path.join(self.experiment_dir, dir, f'{dir}_apo.pdb'))]
		if not compound_dirs:
			raise FileNotFoundError("No compound directories found in the experiment directory.")

		return compound_dirs
		





# Test
if __name__ == '__main__':
	duck = DUckExperiment(
		'/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac_simulations/Experiment',
		wqb_filename='wqb.txt'
		)
	print(duck.compound_dirs)
	duck.summarise_experiment()

	import pandas as pd
	df = pd.read_json('/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac_simulations/Experiment/analysis/experiment_summary.json')

	# Display the DataFrame
	print(df)