"""
Process the output directory containing results of DUck simulations for multiple compounds.
"""



class DUckExperiment():
	"""
	Process the output directory containing results of DUck simulations.
	"""
	def __init__(self, experiment_dir, wqb_filename='wqb.txt'):
		self.experiment_dir = experiment_dir


	def get_simulations(self):
		"""
		Return a list of all simulations in the experiment.
		"""
		pass

	def summarise_experiment(self):
		"""
		Return a json/csv summary of the experiment.
		"""
		pass




# Test
if __name__ == '__main__':
	duck = DUckExperiment('/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac_simulations/Experiment', 'wqb.txt')
	duck.summarise_experiment()
	duck.get_simulations()