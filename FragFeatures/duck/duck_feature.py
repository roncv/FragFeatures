"""
Process the duck output for a feature of a compound.
"""
import os
import json

class DUckFeature():
	"""
	Process the FragFeature metadata for a compound's feature.
	"""
	def __init__(self, feature_dir):
		self.feature_dir = feature_dir


	def get_feature_metadata(self):
		"""
		Return the feature metadata for a compound.
		"""
		pass


	# TODO: Move to utils.py??
	def json_to_dict(self, filename=None):
		"""
		Read in a json file as a dictionary.
		"""
		json_file = os.path.join(self.feature_dir, filename)
		if not os.path.isfile(json_file):
			raise FileNotFoundError(f"`{filename}` not found in the feature directory.")

		# Read the json file
		with open(json_file, 'r') as file:
			data = json.load(file)
		# print(data)
		return data


	def from_pickle(self):
		"""
		Load the features' original `Pose()` object from a pickle file.
		"""
		pass



# Test
if __name__ == '__main__':
	duck_feature = DUckFeature('/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac_simulations/Experiment/cx0270a/A_ARG_144_NH1')
	duck_feature.json_to_dict(filename='protein_feature_metadata.json')
	# duck_feature.get_feature()
	# duck_feature.get_feature_dict()
	# duck_feature.get_feature_json()