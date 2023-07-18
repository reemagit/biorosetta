from setuptools import setup
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README_pypi.md").read_text()

setup(
	name='biorosetta',
	version='0.3.1',
	#packages=find_packages(include=['classes','queries','utils']),
	packages=['biorosetta'],
	url='https://github.com/reemagit/biorosetta',
	author='Enrico Maiorino',
	author_email='enrico.maiorino@gmail.com',
	description='A package to convert gene identifiers between different naming conventions',
	install_requires=['biothings_client','tqdm','pandas'],
	long_description = long_description,
	long_description_content_type = 'text/markdown',
	#package_data={'': ['README_pypi.md']},
	#include_package_data=True
)



