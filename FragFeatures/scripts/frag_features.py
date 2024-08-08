"""
Console Script for FragFeatures
"""
import argparse
# import yaml # NOTE: PyYAML is a dependency??

def hello(hello_who=None):
	'''Hello World function
	'''
	if hello_who is not None:
		print(f"Hello, {hello_who}!")
	else:
		print("Hello, World!")


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
	return args


def parse_input():
	'''Main FragFeatures parser, subparsers define action modes
	'''
	parser = argparse.ArgumentParser(description='Open-source toolkit for extracting fragment features.')
	parser.set_defaults(mode=None)
	modes = parser.add_subparsers(title='Open-source fragment feature extraction toolkit. Choose one of the following actions:', help='', metavar='')

	#Arguments for OpenMM form equilibrated system
	hello = modes.add_parser('hello', help='Hello function from FragFeatures (testing).', description='A simple "Hello World" function.')
	hello.add_argument('-w', '--who', type=str, default = None, help='Say hello to this person.')
	hello.set_defaults(mode='hello')


	args = args_sanitation(parser, modes)

	return args, parser


def main():
	# Parse and sanitize the inputs
	args, parser = parse_input()

	# Chose and perform the specified action
	if args.mode == 'hello':
		hello(hello_who=args.who)
	else:
		parser.print_help()


if __name__ == '__main__':
	main()