import requests
from tqdm import tqdm
import collections
from . import queries
import pandas as pd
from pathlib import Path

def download(url: str, fname: str):
	Path(fname).parent.mkdir(parents=True, exist_ok=True)
	resp = requests.get(url, stream=True)
	total = int(resp.headers.get('content-length', 0))
	# Can also replace 'file' with a io.BytesIO object
	with open(fname, 'wb') as file, tqdm(
		desc=fname,
		total=total,
		unit='iB',
		unit_scale=True,
		unit_divisor=1024,
	) as bar:
		for data in resp.iter_content(chunk_size=1024):
			size = file.write(data)
			bar.update(size)


def make_list(query):
	if (not isinstance(query, collections.abc.Iterable) or isinstance(query, str)):
		return [query]
	return query

def is_list(query):
	return isinstance(query, collections.abc.Iterable) and not isinstance(query, str)

def selection_func(id_out):
	if id_out == 'entr':
		return lambda id_in, id_out_list: min(map(int, id_out_list))
	elif id_out == 'symb':
		return lambda id_in, id_out_list: min(id_out_list, key=len)
	elif id_out == 'ensg':
		return lambda id_in, id_out_list: min(id_out_list)
	else:
		return lambda id_in, id_out_list: min(id_out_list)

def no_intersection(list_of_lists):
	if len(list_of_lists) > 0:
		return len(set.intersection(*[set(elem.split('|')) for elem in list_of_lists])) == 0
	return False

def consensus_elem(list_of_lists):
	flat_list = sum(list_of_lists, [])
	set_list = set(flat_list)
	counts = {elem:flat_list.count(elem) for elem in set_list}
	priorities = {elem: -min([i for i in range(len(list_of_lists)) if elem in list_of_lists[i]]) for elem in set_list}
	orders = {elem: -min([lst.index(elem) for lst in list_of_lists if elem in lst]) for elem in set_list}
	return max(set(flat_list), key=lambda x: (counts[x],priorities[x],orders[x]))

def download_ensembl(path):
	print(f'Downloading to {path}...')
	download(queries.ENSEMBL, path)
	data = pd.read_table(path, header=None, names=['ensg', 'ensp', 'symb', 'synonym', 'entr', 'hgnc'], dtype={'entr': 'str'})[['ensg','ensp','entr','hgnc','symb', 'synonym']]
	data = data[~data.duplicated()]
	data.to_csv(path,sep='\t',index=False)

def download_hgnc(path):
	print(f'Downloading to {path}...')
	download(queries.HGNC, path)
	data = pd.read_table(path, header=0, names=['hgnc', 'symb', 'entr', 'ensg', 'synonym1', 'synonym2'], dtype={'entr': 'str'})[['ensg','entr','hgnc','symb','synonym1', 'synonym2']]
	data = data[~data.duplicated()]
	data.to_csv(path,sep='\t',index=False)