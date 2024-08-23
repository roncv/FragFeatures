"""
Process the output of a compound's DUck simulation.
"""




class DUckCompound():
	"""
	Process the output of a compound's DUck simulation.
	"""
	def __init__(self, experiment_dir, wqb_filename='wqb.txt'):
		self.experiment_dir = experiment_dir


	def validate_simulation(self):
		"""
		Validate that the simulation has run.
		"""
		pass


	# TODO: Make into a property?
	def get_wqb(self):
		"""
		Return the wqb value of a DUck simulation.
		"""
		pass


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



# Test
if __name__ == '__main__':
	duck = DUckCompound('/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac_simulations/Experiment', 'wqb.txt')