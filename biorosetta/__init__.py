from .classes import IDMapper, HGNCBiomartMapper, EnsemblBiomartMapper, MyGeneMapper

id_type_labels = {'ensg':'Ensembl gene ID',
          'entr':'NCBI gene ID (entrezgene)',
          'symb':'Gene symbol',
          'ensp':'Ensembl protein ID',
          'hgnc':'HGNC ID'}

sources = {'ensembl_biomart':
                {'name':'Ensembl Biomart',
                 'url':'http://useast.ensembl.org/biomart/martview',
                 'id_in':['ensg','entr','symb','ensp','hgnc'],
                 'id_out':['ensg','entr','symb','ensp','hgnc']
                },
            'hgnc_biomart':
                {'name':'HGNC Biomart',
                 'url':'https://biomart.genenames.org/',
                 'id_in':['ensg','entr','symb','hgnc'],
                 'id_out':['ensg','entr','symb','hgnc']
                },
            'mygene':
                {'name':'MyGene',
                 'url':'https://mygene.info/',
                 'id_in':['ensg','entr'],
                 'id_out':['ensg','entr','symb']
                }
            }

def list_sources():
    print('ID types:')
    for id_type in id_type_labels:
        print(f'\'{id_type}\' = {id_type_labels[id_type]}')
    print('Sources:')
    for src in sources:
        print(f'\"{src}\": {sources[src]["name"]} ({sources[src]["url"]}) \n\t ID in: {sources[src]["id_in"]} \n\t ID out: {sources[src]["id_out"]}')
    

