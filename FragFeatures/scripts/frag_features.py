"""
Console Script for FragFeatures
"""
import argparse
# import yaml # NOTE: PyYAML is a dependency??
from FragFeatures.utils import timeit

def hello(hello_who=None):
	'''Hello World function
	'''
	if hello_who is not None:
		print(f"Hello, {hello_who}!")
	else:
		print("Hello, World!")


def prepare_duck_experiment(compound_selection, experiment_name, target_dir):
	'''Prepare an experiment for DUck.
	
	Parameters
	----------
	compound_selection : list | str
		Compound selection for the experiment.
	experiment_name : str
		Name of the experiment.
	target_dir : str
		Path to the Fragalysis target's directory.
	'''
	from FragFeatures.duck_input import DUckInput
	# Add some timing
	import time
	start_time = time.time()

	# Prepare the input for DUck simulation
	duck_input = DUckInput(compound_selection=compound_selection,
						   experiment_name=experiment_name,
						   target_dir=target_dir,
						   )
	# print(duck_input.compound_codes)
	duck_input.prepare_experiment()

	# Timing
	end_time = time.time()
	print(f"`prepare_duck_experiment` executed in {end_time - start_time:.4f} seconds")





def args_sanitation(parser, modes):
	'''Sanitize the parser to allow yaml or command line inputs with a proper formating for the rest of the script to work
	'''
	args = parser.parse_args()

	# HELLO
	# check if everything is ok
	if args.mode == 'hello':
		if (args.who is None):
			# This overwrites the function if condition is not met
			modes.choices['hello'].error("You didn't specify who to say hello to.")
			# print("You didn't specify who to say hello to.")
		else:
			pass

	# PREPARE-DUCK
	# TODO: Type check for compound_selection
	elif args.mode == 'prepare-duck':
		print(args.compound_selection, args.experiment_name, args.target_dir)
		if (args.compound_selection is None) or (args.experiment_name is None) or (args.target_dir is None):
			# This overwrites the function if condition is not met
			modes.choices['prepare-duck'].error("You didn't specify all the required arguments.")
			# print("You didn't specify all the required arguments.")
	else:
		pass

	return args


def parse_input():
	'''Main FragFeatures parser, subparsers define action modes
	'''
	parser = argparse.ArgumentParser(description='Open-source toolkit for extracting fragment features.')
	parser.set_defaults(mode=None)
	modes = parser.add_subparsers(title='Open-source fragment feature extraction toolkit. Choose one of the following actions:', help='', metavar='')


	#Arguments for hello function (hello)
	hello = modes.add_parser('hello', help='Hello function from FragFeatures (testing).', description='A simple "Hello World" function.')
	hello.add_argument('-w', '--who', type=str, default = None, help='Say hello to this person.')
	hello.set_defaults(mode='hello')


	#Arguments for preparing a DUck experiment (prepare-duck)
	prepare_duck = modes.add_parser('prepare-duck', help='Prepare an experiment for DUck.')
	prepare_duck.add_argument('-c', '--compound-selection', type=str, nargs='+', default = None, help='Compound selection for the experiment.')
	prepare_duck.add_argument('-e', '--experiment-name', type=str, default = None, help='Name of the experiment.')
	prepare_duck.add_argument('-t', '--target-dir', type=str, default = None, help='Path to the Fragalysis target\'s directory.')
	# prepare_duck.add_argument('-a', '--all-compounds', action='store_true', help='Use all compounds in the target.')
	prepare_duck.set_defaults(mode='prepare-duck')

	args = args_sanitation(parser, modes)
	# print(args)

	return args, parser



@timeit
def main():
	# Parse and sanitize the inputs
	args, parser = parse_input()

	# Chose and perform the specified action
	if args.mode == 'hello':
		hello(hello_who=args.who)
	elif args.mode == 'prepare-duck':
		prepare_duck_experiment(compound_selection=args.compound_selection,
								experiment_name=args.experiment_name,
								target_dir=args.target_dir)
	else:
		parser.print_help()


if __name__ == '__main__':
	main()