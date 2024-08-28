"""
Various utility functions for FragFeatures.
"""


# Decorator to time functions
def timeit(func):
	"""
	Decorator to time functions.
	"""
	import time

	def wrapper(*args, **kwargs):
		start = time.time()
		result = func(*args, **kwargs)
		end = time.time()
		print(f"\nfunction {func.__name__} took {end - start:.2f}s to run.")
		return result
	return wrapper

def dict_to_json(data, filename, mode='w', indent=4):
	"""
	Write a dictionary to a JSON file.
	"""
	import json
	with open(filename, f'{mode}') as f:
		json.dump(data, f, indent=indent)