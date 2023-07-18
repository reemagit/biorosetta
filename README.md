# BioRosetta

This is a package to map gene names between different naming conventions.

Motivation: while there are popular packages for gene identifier mapping (e.g. ENSG, NCBI, HGNC) in the R environment (e.g. AnnotationDB), there is no standard solution in python.

## Package Features

```python
import biorosetta as br
```

- **Source-based system**: Instead of relying on a single repository for mapping gene identifiers, biorosetta integrates results from different repositories, or "sources". Biorosetta supports two types of sources under the same interface. (1) *Local sources*: biorosetta downloads a local version of the conversion tables from popular repositories (Ensembl Biomart, HGNC Biomart). Best option for highly reproducible gene conversion outputs that do not change over time (e.g. for scientific article preparation). (2) *Remote sources*: biorosetta interfaces to remote web service applications (MyGene) to convert gene names. Best option for highly up-to-date conversion.
- **Priority system**: The user can specify an order of source priority, so that when different sources produce a different conversion output it is possible to define the most trusted result. 
- **Conversion report**: For critical gene mapping applications, biorosetta can optionally generate a report table that specifies the mapping output of each separate source and highlight where there have been mismatches between outputs, so that these mapping results can be investigated further.
- **Multi-hits policies**: When multiple possible mapping outputs ("hits") are found, one can choose the policy for integrating them: "all": concatenates all the ID outputs with a pipe ("|") symbol (e.g. "foo|bar|baz"). "consensus": compares the output hits across different sources to select the ID that appears most frequently across different sources. E.g. if source A outputs "foo|bar" and source B outputs "bar|baz" then "bar" is selected as the final output.
- **Gene name synonyms**: Gene symbols have many synonyms. When available within the source, biorosetta integrates the synonym information for one-way mapping from gene symbol to any other ID.
- **Speed**: Biorosetta performs vectorized operations to achieve maximum efficiency. 

## Usage

```python
import biorosetta as br
```

List current supported sources and gene ID types:
```python
br.list_sources()
```
Function output:
```
ID types:
'ensg' = Ensembl gene ID
'entr' = NCBI gene ID (entrezgene)
'symb' = Gene symbol
'ensp' = Ensembl protein ID
'hgnc' = HGNC ID
Sources:
"ensembl_biomart": Ensembl Biomart (http://useast.ensembl.org/biomart/martview) 
	 ID in: ['ensg', 'entr', 'symb', 'ensp', 'hgnc'] 
	 ID out: ['ensg', 'entr', 'symb', 'ensp', 'hgnc']
"hgnc_biomart": HGNC Biomart (https://biomart.genenames.org/) 
	 ID in: ['ensg', 'entr', 'symb', 'hgnc'] 
	 ID out: ['ensg', 'entr', 'symb', 'hgnc']
"mygene": MyGene (https://mygene.info/) 
	 ID in: ['ensg', 'entr'] 
	 ID out: ['ensg', 'entr', 'symb']
```

## Basic usage

This is the basic command that uses all the sources (local and remote) ordered with the following order of priority:
1. Ensembl Biomart
2. HGNC Biomart
3. MyGene

```python
idmap = br.IDMapper('all') # Equivalent to br.IDMapper([br.EnsemblBiomartMapper(),br.HGNCBiomartMapper(),br.MyGeneMapper()])
idmap.convert('ENSG00000159388','ensg','entr') # Outputs '7832'
```

Note that if a list of gene IDs is provided, the function returns a list:

```python
idmap.convert(['ENSG00000159388','ENSG00000121022','ENSG00000134574'],'ensg','entr') # Outputs ['7832', '10987', '1643']
```

## Examples

### Single source mapping

Here we map ENSG to NCBI ID (entrezgene) using only Ensembl Biomart as reference data:
```python
idmap = br.IDMapper([br.EnsemblBiomartMapper()]) # Single local source
idmap.convert('ENSG00000159388','ensg','entr')
```

### Multiple sources mapping 

We map ENSG to NCBI ID (entrezgene) using only Ensembl Biomart as reference data:
```python
idmap = br.IDMapper('all') # Equivalent to br.IDMapper([br.EnsemblBiomartMapper(),br.HGNCBiomartMapper(),br.MyGeneMapper()])
idmap.convert('ENSG00000271254','ensg','entr') # Returns '102724250'
```
However, HGNC Biomart does not have the conversion. So if we mapped only with HGNC Biomart we would get

```python
idmap = br.IDMapper('hgnc_biomart') # Equivalent to br.IDMapper([br.HGNCBiomartMapper()])
idmap.convert('ENSG00000271254','ensg','entr') # Returns 'N/A'
```

__Note: when using **multi_hits = "first"** in the **convert()** function, the order of the sources in the list provided to IDMapper defines the priority of each source in the mapping output. The output of each source will be the first matching ID returned, and the final conversion output will be the output of the source with the highest priority that found a match for the input ID. See section "Multi hits" for more information.__

### Generation of report of conversion results

For each conversion, we can choose to generate a report of all the conversion results across each source, to be able to check single responses and diagnose bad conversions. Usually convenient for conversions of a small number of gene IDs in situations where accuracy is critical. 

```python
idmap = br.IDMapper([br.EnsemblBiomartMapper(),br.HGNCBiomartMapper(),br.MyGeneMapper()])
idmap.convert(gene_list[:3],'entr','ensg',df=True)
```
These instructions return the following pandas DataFrame:

|    |   input | output          | ensembl         | hgnc            | mygene          | mismatch   |   ensembl_hits |   hgnc_hits |   mygene_hits |
|---:|--------:|:----------------|:----------------|:----------------|:----------------|:-----------|---------------:|------------:|--------------:|
|  0 |    1643 | ENSG00000134574 | ENSG00000134574 | ENSG00000134574 | ENSG00000134574 | False      |              1 |           1 |             1 |
|  1 |   10987 | ENSG00000121022 | ENSG00000121022 | ENSG00000121022 | ENSG00000121022 | False      |              1 |           1 |             1 |
|  2 |   23369 | ENSG00000055917 | ENSG00000055917 | ENSG00000055917 | ENSG00000055917 | False      |              1 |           1 |             1 |

Columns:
- 'input': the queried gene ID list
- 'output': the final ID conversion output (returned also if df=False)
- 'ensembl','hgnc', 'mygene': report the conversion for each source
- 'mismatch': whether there is a mismatch between sources conversion (different IDs returned)
- '{source}_hits': number of possible output ID conversions found by each source.

### Multi hits

There are two types of policies we have to define to convert an input gene ID with multiple sources:
- Each source can return more than one output ID for each input ID, so we have to choose what ID will be return among the list
- The final outputs of all the sources have to be integrated to return a single integrated IDMapper output

These policies are handled through the "multi_hits" argument of the convert() function of IDMapper. Here we define the policies available:

##### 1. Multi-hits = 'first'

For each source, only the first conversion hit will be the source output. The outputs of different sources are integrated by selecting the output of the source with (1) the highest priority and (2) a successful conversion (i.e. not 'N/A').
```python
idmap = br.IDMapper([br.MyGeneMapper(),br.EnsemblBiomartMapper(),br.HGNCBiomartMapper()], multi_hits='first', df=True) # Notice that MyGene is first
idmap.convert(['ENSG00000210049','ENSG00000211459','ENSG00000210082'],'ensg','symb',df=True) # Outputs ['TRNF', 'RNR1', 'RNR2']
```

And this is the final conversion report. Note that the MyGene source disagrees with the Ensembl Biomart and HGNC Biomart, but it has highest priority.

|    | input           | output   | mygene   | ensembl   | hgnc    | mismatch   |   mygene_hits |   ensembl_hits |   hgnc_hits |
|---:|:----------------|:---------|:---------|:----------|:--------|:-----------|--------------:|---------------:|------------:|
|  0 | ENSG00000210049 | TRNF     | TRNF     | MT-TF     | MT-TF   | True       |             1 |              1 |           1 |
|  1 | ENSG00000211459 | RNR1     | RNR1     | MT-RNR1   | MT-RNR1 | True       |             1 |              1 |           1 |
|  2 | ENSG00000210082 | RNR2     | RNR2     | MT-RNR2   | MT-RNR2 | True       |             1 |              1 |           1 |

##### 2. multi_hits = 'consensus'

Consensus mapping returns the most frequent ID returned across all the sources. 
```python
idmap = br.IDMapper([br.MyGeneMapper(),br.EnsemblBiomartMapper(),br.HGNCBiomartMapper()]) # Multiple sources
idmap.convert(['ENSG00000210049','ENSG00000211459','ENSG00000210082'],'ensg','symb',df=True, multi_hits='consensus') # Returns ['MT-TF', 'MT-RNR1', 'MT-RNR2']
```

Note that the final output are not the IDs returned by MyGene in our previous example since they were not the most frequent outputs. Here is the conversion report.

|    | input           | output   | mygene   | ensembl   | hgnc    | mismatch   |   mygene_hits |   ensembl_hits |   hgnc_hits |
|---:|:----------------|:---------|:---------|:----------|:--------|:-----------|--------------:|---------------:|------------:|
|  0 | ENSG00000210049 | MT-TF    | TRNF     | MT-TF     | MT-TF   | True       |             1 |              1 |           1 |
|  1 | ENSG00000211459 | MT-RNR1  | RNR1     | MT-RNR1   | MT-RNR1 | True       |             1 |              1 |           1 |
|  2 | ENSG00000210082 | MT-RNR2  | RNR2     | MT-RNR2   | MT-RNR2 | True       |             1 |              1 |           1 |



##### 2. multi_hits = 'all'

Return all ID outputs for each source concatenated with a pipe ("|") symbol. 

```python
idmap = br.IDMapper([br.EnsemblBiomartMapper(),br.HGNCBiomartMapper(),br.MyGeneMapper()]) # Multiple sources
idmap.convert(['ENSG00000159388','ENSG00000211459','ENSG00000159388'],'ensg','ensp', multi_hits='all') # Returns ['ENSP00000433553|ENSP00000290551', 'N/A', 'ENSP00000433553|ENSP00000290551']
```

### Setting default value for failed conversions (e.g. "N/A")

The default value for failed ID conversions it "N/A". This value can be modified when creating the IDMapper class.

```python
idmap = br.IDMapper('all')
idmap.convert(['imaginary_gene1','imaginary_gene2'],'ensg','entr') # Returns ['N/A', 'N/A']
```

If we set "fill_value='passthrough'" the input ID will be returned when the output ID is not found:

```python
idmap = br.IDMapper('all', fill_value='passthrough')
idmap.convert(['imaginary_gene1','imaginary_gene2'],'ensg','entr') # Returns ['imaginary_gene1', 'imaginary_gene2']
```


## Currently implemented sources and gene identifiers

__Note: Not all gene identifiers work on all sources, depending on availability__

#### Sources:
- Ensembl Biomart (https://useast.ensembl.org/biomart/martview): Local source
- HGNC Biomart (GeneNames) (https://biomart.genenames.org/): Local source
- MyGene (https://mygene.info): Remote source

#### Gene Identifiers:

- "ensg": Ensembl gene ID (all sources)
- "entr": NCBI gene ID (entrezgene, all sources)
- "symb": Gene name (symbol, all sources)
- "ensp": Ensembl protein ID (ENSP, Ensembl Biomart only)
- "hgnc": HGNC ID (Ensembl Biomart and HGNC Biomart only)



## IDMapper initialization aliases

```python
idmap = br.IDMapper('all') # equivalent to [br.EnsemblBiomartMapper(),br.HGNCBiomartMapper(),br.MyGeneMapper()]
idmap = br.IDMapper('local') # equivalent to [br.EnsemblBiomartMapper(),br.HGNCBiomartMapper()]
idmap = br.IDMapper('remote') # equivalent to [br.MyGeneMapper()]
idmap = br.IDMapper('ensembl_biomart') # equivalent to [br.EnsemblBiomartMapper()]
idmap = br.IDMapper('hgnc_biomart') # equivalent to [br.HGNCBiomartMapper()]
idmap = br.IDMapper('mygene') # equivalent to [br.MyGene()]
```

# Acknowledgments

Thanks to David Deritei for the help with design and debugging.