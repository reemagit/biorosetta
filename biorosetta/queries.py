ENSEMBL='''http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_peptide_id" /><Attribute name = "external_gene_name" /><Attribute name = "external_synonym" /><Attribute name = "entrezgene_id" /><Attribute name = "hgnc_id" /></Dataset></Query>'''
HGNC='''http://biomart.genenames.org/martservice/results?query=<!DOCTYPE Query><Query client="biomartclient" processor="TSV" limit="-1" header="1"><Dataset name="hgnc_gene_mart" config="hgnc_gene_config"><Filter name="hgnc_gene__status_1010" value="Approved" filter_list=""/><Attribute name="hgnc_gene__hgnc_gene_id_1010"/><Attribute name="hgnc_gene__approved_symbol_1010"/><Attribute name="hgnc_gene__ncbi_gene__gene_id_1026"/><Attribute name="hgnc_gene__ensembl_gene__ensembl_gene_id_104"/><Attribute name="hgnc_gene__hgnc_alias_symbol__alias_symbol_108"/><Attribute name="hgnc_gene__hgnc_previous_symbol__previous_symbol_1012"/></Dataset></Query>'''