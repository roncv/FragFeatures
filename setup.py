from setuptools import setup, find_packages
from pathlib import Path


# # Versioning
# def _get_version() -> str:
# 	"""Read FragFeatures/VERSION.txt and return its contents."""
# 	path = Path("PyCEC").resolve()
# 	version_file = path / "VERSION.txt"
# 	return version_file.read_text().strip()

# TODO: Add versioning
# TODO: Transition to pyproject.toml for versioningß

# version = _get_version()

requirements = [
	# 'numpy',
	# 'MDAnalysis',
	# 'tqdm',
	# 'matplotlib'
    'rdkit',
    # 'molparse',
    # 'pandas',
    'openmm',
    'prolif',
	]

# TODO: Add long description linked to README.md

setup(
	name='FragFeatures',
	version=1.0,
	author='Ronald Cvek',
	author_email='ronaldcvek@gmail.com',
	description=(
		'Tool to extract protein-ligand interaction '
		'fingerprints from Fragalysis data.'
		),
	url='https://github.com/roncv/FragFeatures',
	packages=find_packages(),
	install_requires=requirements,
	classifiers=[
		'Development Status :: 3 - Alpha',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 3',
		'Topic :: Scientific/Engineering :: Chemistry',
	],
	python_requires='>=3.10',
	# package_data={  # Optional: include data files
	#     'your_package_name': ['data/*']
	# },
	entry_points={  # Optional: define command-line scripts
		'console_scripts': [
			'fragfeat=FragFeatures.scripts.frag_features:main',
		],
	},
)
