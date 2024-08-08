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
			  compound_code: list,
			  experiment_name: str,
			  target_dir: str,
			  ):
		self.compound_code = compound_code
		self.experiment_name = experiment_name
		self.target_dir = target_dir
		self.target = TargetParser(target_dir)


	def get_compound(self, code):
		pass










if __name__ == '__main__':
	# Testing
	print('Testing...')
	# pass