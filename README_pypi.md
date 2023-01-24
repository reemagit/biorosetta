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

## Currently implemented sources and gene identifiers

Sources:
- Ensembl Biomart (https://useast.ensembl.org/biomart/martview): Local source
- HGNC Biomart (GeneNames) (https://biomart.genenames.org/): Local source
- MyGene (https://mygene.info): Remote source

Gene Identifiers:

- "ensg": Ensembl gene ID (all sources)
- "entr": NCBI gene ID (entrezgene, all sources)
- "symb": Gene name (symbol, all sources)
- "ensp": Ensembl protein ID (ENSP, Ensembl Biomart only)
- "hgnc": HGNC ID (Ensembl Biomart and HGNC Biomart only)

## Usage

See up-to-date documentation and examples at [github repo](https://www.google.com).

