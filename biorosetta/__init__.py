from .classes import IDMapper, HGNCBiomartMapper, EnsemblBiomartMapper, MyGeneMapper


def get_mapper(source='all'):
	if source == 'ensembl_biomart':
		return IDMapper([EnsemblBiomartMapper()])
	elif source == 'hgnc_biomart':
		return IDMapper([HGNCBiomartMapper()])
	elif source == 'mygene':
		return IDMapper([MyGeneMapper()])
	elif source == 'all':
		return IDMapper([EnsemblBiomartMapper(),HGNCBiomartMapper(),MyGeneMapper()])
	elif source == 'local':
		return IDMapper([EnsemblBiomartMapper(),HGNCBiomartMapper()])
	elif source == 'remote':
		return IDMapper([MyGeneMapper()])
	else:
		raise ValueError(f'Source specified ({source}) is invalid')