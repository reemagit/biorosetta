from biothings_client import get_client
from pathlib import Path
from .utils import *
import os
import pickle
import numpy as np

lib_folder = os.path.dirname(os.path.realpath(__file__))



class Source:
	def __init__(self, source_id, fill_value='N/A'):
		self.source_id = source_id
		self._fill_value = fill_value

	def sanitize(self, id_list, id_in, id_out):
		if not self.has_id_in_type(id_in):
			raise ValueError(f'{self.source_id}: invalid input ID type')
		if not self.has_id_out_type(id_out):
			raise ValueError(f'{self.source_id}: invalid output ID type')

		return list(map(str, make_list(id_list)))

	def filter_multi_hits(self, out_df, multi_hits):
		if multi_hits == 'first':
			out_df[out_df.str.contains('\|')] = out_df[out_df.str.contains('\|')].str.split('|').str[0]
		if multi_hits == 'shortest':
			out_df[out_df.str.contains('\|')] = out_df[out_df.str.contains('\|')].str.split('|').apply(lambda x: min(x,key=len))
		elif callable(multi_hits):
			out_df[out_df.str.contains('\|')] = map(multi_hits, out_df[out_df.str.contains('\|')].str.split('|').tolist())
		elif multi_hits == 'all': # just here as placeholder, but any other value will do
			pass
		return out_df


class LocalSource(Source):
	def __init__(self, source_id, data, fill_value='N/A'):
		super().__init__(source_id, fill_value=fill_value)

		self._data = data
		self.init()


	def init(self):
		data = self._data
		self._lookup = {}

		id_types = data.columns.tolist()
		for id_in in id_types:
			self._lookup[id_in] = {}
			for id_out in id_types:
				if id_in == id_out:
					continue

				self._lookup[id_in][id_out] = self.gen_lookup_table(data, id_in, id_out,
																	lambda id_in, id_out_list: '|'.join(
																		id_out_list))


	def convert(self, id_list, id_in, id_out, multi_hits='first', df=False):
		'''
			:param list id_list: List of IDs to map
			:param str id_in: Input ID type
			:param str id_out: Output ID type
			:param str multi_hits: Aggregation method for multi hits. Possible values:
				- 'first': Returns the first ID returned by each source, and selects the ID from the source with highest priority
				- 'consensus': Returns the most occurring ID returned by all the sources
				- 'all': Returns all IDs returned by each source, separated by a pipe ('|') symbol
			:param bool df: Whether to return a DataFrame with a full reports of the sources responses (True) or just the converted IDs (False)
			:return: List of IDs if df==True, or DataFrame otherwise
			:rtype: list or DataFrame
		'''
		multi_ids = is_list(id_list)
		id_list = self.sanitize(id_list, id_in, id_out)
		out_df = self._lookup[id_in][id_out].reindex(id_list, fill_value=self._fill_value)
		out_df = self.filter_multi_hits(out_df, multi_hits)
		if self._fill_value == 'passthrough':
			out_df = pd.Series(np.where(out_df == 'passthrough', out_df.index, out_df.values), index=out_df.index)
		if df:
			return out_df
		else:
			output = out_df.values
			if multi_ids:
				return output.tolist()
			else:
				return output.tolist()[0]


	def gen_lookup_table(self, data, id_in, id_out, selection_func):
		pruned = data[data[id_in].notna() & data[id_out].notna()]
		pruned = pruned[~pruned.duplicated(subset=[id_in, id_out])]
		id_in_unq = ~pruned.duplicated(subset=[id_in], keep=False)
		unq = pruned.loc[id_in_unq].set_index(id_in)[id_out]

		if (~id_in_unq).sum() > 0:
			nonunq = pruned.loc[~id_in_unq, [id_in, id_out]].groupby(id_in).apply(
				lambda x: selection_func(x.name, x[id_out]))
			# nonunq = pruned.loc[~id_in_unq, [id_in, id_out]].groupby(id_in).apply(lambda x: print(x))
			lookup = pd.concat([unq, nonunq])
		else:
			lookup = unq
		return lookup

	def integrate_synonyms(self, data, id_orig, id_synonym):
		pruned = data[data[id_orig].notna() & data[id_synonym].notna()]
		pruned = pruned[~pruned.duplicated(subset=[id_orig, id_synonym])]
		pruned = pruned[~pruned.duplicated(subset=[id_synonym], keep=False)]  # we just remove inconsistent synonyms
		for id_out in self._lookup[id_orig].keys():
			lookup = self._lookup[id_orig][id_out]
			pruned_curr = pruned.loc[~pruned[id_synonym].isin(lookup.index) & pruned[id_orig].isin(lookup.index)].copy()
			pruned_curr[id_out] = lookup.loc[pruned_curr[id_orig].values].values
			if pruned_curr.shape[0] > 0:
				pruned_curr = pruned_curr.set_index(id_synonym)[id_out].squeeze()
				self._lookup[id_orig][id_out] = pd.concat([lookup, pruned_curr])

	def has_id_in_type(self, id_type):
		return id_type in self._data.columns

	def has_id_out_type(self, id_type):
		return id_type in self._data.columns

	def lookup_size(self):
		for id_in, id_in_dict in self._lookup.items():
			for id_out, id_out_dict in id_in_dict.items():
				print(id_in, ' ->', id_out, ':', len(id_out_dict))

	def build_cache(self, cache_path):
		with open(cache_path, 'wb') as f:
			pickle.dump([self._data, self._lookup], f)

	def load_cache(self, cache_path):
		with open(cache_path, 'rb') as f:
			data,lookup = pickle.load(f)
		self._data = data
		self._lookup = lookup


class EnsemblBiomartMapper(LocalSource):
	def __init__(self, data_path=None, symb_aliases=True, fill_value='N/A'):
		'''
			:param str data_path: Path where to store the local data for this source (default: internal package folder)
			:param bool symb_aliases: Whether to download and integrate the symbol aliases and synonyms in the dictionary
			:param str fill_value: Default value returned when the ID is not found in the source
		'''
		self._source_label = 'Ensembl Biomart'
		if data_path is None:
			data_path = lib_folder + '/data/ensembl.tsv'
		cache_path = data_path.replace('.tsv', '.pickle')
		if Path(cache_path).exists():
			print('- Loading lookup tables from cache (use function EnsemblBiomartMapper.download_data() to force new download)')
			self.load_cache(cache_path)
			self.source_id = 'ensembl'
			self._fill_value = fill_value
		else:
			if not Path(data_path).exists():
				print('- Biomart data has not been downloaded yet.')
				print(f'- Downloading up-to-date gene annotation data from Ensembl Biomart (http://www.ensembl.org/biomart) to {data_path}.')
				print('- This operation has to be performed only once and it lasts less than few minutes.')
				EnsemblBiomartMapper.download_data(data_path)
				print('- Download completed')
			data = pd.read_table(data_path, sep='\t', dtype={'entr': 'str'})
			main_data = data[['ensg', 'ensp', 'entr', 'hgnc', 'symb']]
			super().__init__(source_id='ensembl', data=main_data, fill_value=fill_value)
			if symb_aliases:
				syn_data = data[['symb', 'synonym']]
				self.integrate_synonyms(syn_data, 'symb', 'synonym')
			self.build_cache(cache_path)

	@staticmethod
	def download_data(data_path=None):
		if data_path is None:
			data_path = lib_folder + '/data/ensembl.tsv'
		download_ensembl(data_path)



class HGNCBiomartMapper(LocalSource):
	def __init__(self, data_path=None, symb_aliases=True, fill_value='N/A'):
		'''
			:param str data_path: Path where to store the local data for this source (default: internal package folder)
			:param bool symb_aliases: Whether to download and integrate the symbol aliases and synonyms in the dictionary
			:param str fill_value: Default value returned when the ID is not found in the source
		'''
		self._source_label = 'HGNC Biomart'
		if data_path is None:
			data_path = lib_folder + '/data/hgnc.tsv'
		cache_path = data_path.replace('.tsv', '.pickle')
		if Path(cache_path).exists():
			print('- Loading lookup tables from cache (use function HGNCBiomartMapper.download_data() to force new download)')
			self.load_cache(cache_path)
			self.source_id = 'hgnc'
			self._fill_value = fill_value
		else:
			if not Path(data_path).exists():
				print('- Biomart data has not been downloaded yet.')
				print(f'- Downloading up-to-date gene annotation data from HGNC Biomart (http://biomart.genenames.org) to {data_path}.')
				print('- This operation has to be performed only once and it lasts less than few minutes.')
				HGNCBiomartMapper.download_data(data_path)
				print('- Download completed')
			data = pd.read_table(data_path, sep='\t', dtype={'entr': 'str'})
			main_data = data[['ensg', 'entr', 'hgnc', 'symb']]
			syn_data = data[['symb', 'synonym1', 'synonym2']]
			super().__init__(source_id='hgnc', data=main_data, fill_value=fill_value)
			if symb_aliases:
				self.integrate_synonyms(syn_data, 'symb', 'synonym1')
				self.integrate_synonyms(syn_data, 'symb', 'synonym2')
			self.build_cache(cache_path)

	@staticmethod
	def download_data(data_path=None):
		if data_path is None:
			data_path = lib_folder + '/data/hgnc.tsv'
		download_hgnc(data_path)

class RemoteSource(Source):
	def __init__(self, source_id, fill_value='N/A'):
		super().__init__(source_id, fill_value=fill_value)


class MyGeneMapper(RemoteSource):
	mg = get_client('gene')

	def __init__(self, fill_value='N/A'):
		'''
			:param str fill_value: Default value returned when the ID is not found in the source
		'''
		self._source_label = 'MyGene'
		super().__init__('mygene', fill_value=fill_value)

	def convert(self, id_list, id_in, id_out, multi_hits='first', df=False):
		'''
				:param list id_list: List of IDs to map
				:param str id_in: Input ID type
				:param str id_out: Output ID type
				:param str multi_hits: Aggregation method for multi hits. Possible values:
					- 'first': Returns the first ID returned by each source, and selects the ID from the source with highest priority
					- 'consensus': Returns the most occurring ID returned by all the sources
					- 'all': Returns all IDs returned by each source, separated by a pipe ('|') symbol
				:param bool df: Whether to return a DataFrame with a full reports of the sources responses (True) or just the converted IDs (False)
				:return: List of IDs if df==True, or DataFrame otherwise
				:rtype: list or DataFrame
		'''
		id_list = self.sanitize(id_list, id_in, id_out)

		id_relabel = {'entr': 'entrezgene', 'symb': 'symbol', 'ensg': 'ensembl.gene', 'ensp':'ensembl.protein', 'hgnc':'HGNC'}
		output = MyGeneMapper.mg.getgenes(id_list, scopes=id_relabel[id_in], fields=id_relabel[id_out],
									species='human', as_dataframe=True, returnall=False)

		if output.shape[1]==1 and 'notfound' in output.columns:
			out_df = pd.Series([self._fill_value] * len(id_list), index=id_list)
		else:
			out_df = output[id_relabel[id_out]].fillna(self._fill_value)
			out_df = self.filter_multi_hits(out_df, multi_hits)
			out_df = out_df[~out_df.index.duplicated()].reindex(id_list) # to remove duplicated input IDs in query
		if self._fill_value == 'passthrough':
			out_df = pd.Series(np.where(out_df == 'passthrough', out_df.index, out_df.values), index=out_df.index)
		out_df = out_df.astype(str)
		if df:
			return out_df
		else:
			return out_df.values

	def has_id_in_type(self, id_type):
		return id_type in ['entr', 'ensg']

	def has_id_out_type(self, id_type):
		return id_type in ['entr', 'ensg', 'symb']


class IDMapper:
	def __init__(self, sources, fill_value=None):
		'''
			:param list or str sources: List of source objects or string with possible values:
				- 'ensembl_biomart': Ensembl Biomart source (local)
				- 'hgnc_biomart': HGNC Biomart source (local)
				- 'mygene': MyGene source (remote)
				- 'all': All sources with following priority list: [Ensembl Biomart,HGNC Biomart, MyGene]
				- 'local': All local sources (Ensembl Biomart, HGNC Biomart)
				- 'remote': All remote sources (MyGene)
			:param str fill_value: Value to return when output ID is not found. 
				- 'N/A': (Default)
				- 'passthrough': Return input ID
				- None: Inherit existing fill values from specified sources
				- any other string: Return any other string
			:return: IDMapper object
			:rtype: IDMapper
		'''
		sources = IDMapper.get_sources(sources)
		self._sources = make_list(sources)
		self._src_ids = [src.source_id for src in self._sources]
		if fill_value is None:
			self._fill_value = self._sources[0]._fill_value
		else:
			self._fill_value = fill_value
			for i in range(len(self._sources)):
				self._sources[i]._fill_value = fill_value
		

	@staticmethod
	def get_sources(source='all'):

		if isinstance(source, str):
			if source == 'ensembl_biomart':
				return [EnsemblBiomartMapper()]
			elif source == 'hgnc_biomart':
				return [HGNCBiomartMapper()]
			elif source == 'mygene':
				return [MyGeneMapper()]
			elif source == 'all':
				return [EnsemblBiomartMapper(), HGNCBiomartMapper(), MyGeneMapper()]
			elif source == 'local':
				return [EnsemblBiomartMapper(), HGNCBiomartMapper()]
			elif source == 'remote':
				return [MyGeneMapper()]
			else:
				raise ValueError(f'Source specified ({source}) is invalid')
		else:
			return source

	def get_source(self, source_id):
		return self._sources[self._src_ids.index(source_id)]

	def convert(self, id_list, id_in, id_out, multi_hits='first', df=False):
		'''
			:param list id_list: List of IDs to map
			:param str id_in: Input ID type
			:param str id_out: Output ID type
			:param str multi_hits: Aggregation method for multi hits. Possible values:
				- 'first': Returns the first ID returned by each source, and selects the ID from the source with highest priority
				- 'consensus': Returns the most occurring ID returned by all the sources
				- 'all': Returns all IDs returned by each source, separated by a pipe ('|') symbol
			:param bool df: Whether to return a DataFrame with a full reports of the sources responses (True) or just the converted IDs (False)
			:return: List of IDs if df==True, or DataFrame otherwise
			:rtype: list or DataFrame
		'''
		multi_ids = is_list(id_list)
		id_list = make_list(id_list)
		src_ids = [src_id for src_id, src in zip(self._src_ids, self._sources) if src.has_id_in_type(id_in) and src.has_id_out_type(id_out)]
		if len(src_ids) < len(self._src_ids):
			print('One or more sources do not support the requested input/output ID type mapping: {}'.format(','.join(src_id for src_id in self._src_ids if src_id not in src_ids)))
			print('Mapping will be executed using the following source(s): {}'.format(','.join(src_ids)))
		if len(src_ids) == 0:
			raise ValueError('Input or output ID type not supported by selected sources')
		out_df = pd.DataFrame([], columns=src_ids)
		out_df['input'] = id_list
		for src_id in src_ids:
			out_df[src_id] = self.get_source(src_id).convert(id_list, id_in, id_out, multi_hits='all' if multi_hits == 'consensus' else multi_hits, df=False)

		if multi_hits == 'consensus':
			def get_consensus(x):
				if (x != self._fill_value).sum() > 0:
					return consensus_elem(x[x != self._fill_value].str.split('|').tolist())
				else:
					return self._fill_value

			out_df['output'] = out_df[src_ids].apply(get_consensus,axis=1)
			#id_list_out = out_df['consensus']
		else:
			id_list_out = out_df[src_ids[0]]
			for i in range(1, len(src_ids)):
				idx = id_list_out == self._fill_value
				if idx.sum().squeeze() == 0:
					break
				id_list_out[idx] = out_df.loc[idx, src_ids[i]]
			out_df['output'] = id_list_out
		out_df.columns = out_df.columns.get_level_values(0)
		if df:
			if multi_hits == 'consensus' or multi_hits == 'all':
				out_df['mismatch'] = out_df[src_ids].apply(lambda x: no_intersection(x[x != self._fill_value].values),axis=1)
			else:
				out_df['mismatch'] = out_df[src_ids].apply(lambda x: x[x != self._fill_value].nunique() > 1, axis=1)
			for src_id in src_ids:
				out_df.loc[out_df[src_id] != 'N/A', f'{src_id}_hits'] = out_df[src_id].str.split('|').apply(lambda x: len(x))
				out_df.loc[out_df[src_id] == 'N/A', f'{src_id}_hits'] = 0
			return out_df[['input','output'] + [col for col in out_df.columns.tolist() if col not in ['input','output']]]
		else:
			if multi_ids:
				return out_df['output'].tolist()
			else:
				return out_df['output'].tolist()[0]

	def entr2ensg(self, id_list, df=False):
		return self.convert(id_list, id_in='entr', id_out='ensg', df=df)

	def entr2symb(self, id_list, df=False):
		return self.convert(id_list, id_in='entr', id_out='symb', df=df)

	def ensg2entr(self, id_list, df=False):
		return self.convert(id_list, id_in='ensg', id_out='entr', df=df)

	def ensg2symb(self, id_list, df=False):
		return self.convert(id_list, id_in='ensg', id_out='symb', df=df)

	def symb2entr(self, id_list, df=False):
		return self.convert(id_list, id_in='symb', id_out='entr', df=df)

	def symb2ensg(self, id_list, df=False):
		return self.convert(id_list, id_in='symb', id_out='ensg', df=df)
