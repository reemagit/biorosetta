from setuptools import setup, find_packages

setup(
	name='biorosetta',
	version='0.1.4',
	#packages=find_packages(include=['classes','queries','utils']),
	packages=['biorosetta'],
	url='https://github.com/reemagit/biorosetta',
	author='Enrico Maiorino',
	author_email='enrico.maiorino@gmail.com',
	description='A package to convert gene identifiers between different naming conventions',
    install_requires=['biothings_client','tqdm','pandas']
)


