"""
Main script for parsing through a Fragalysis target's directory and extract compound/fragment fingerprints.

Inspired by 
"""
if __name__ == '__main__':
	# Conditional imports only when running as the main script
	from parse_target import TargetParser, Pose
else:
	from FragFeatures.parse_target import TargetParser, Pose

import molparse as mp

import numpy as np
import pandas as pd

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


	def get_compound(self, compound_code):
		pass


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

	def get_compound(self, compound_code):
		"""
		Return a compound's metadata with the given code.
		"""
		








if __name__ == '__main__':
	# Testing
	print('Testing...')
	duck_input = DUckInput(compound_selection=['cx0270a', 'cx0281a'], experiment_name='test', target_dir='/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac')
	# duck_input.validate_compounds()
	print(duck_input.compound_codes)
	# pass