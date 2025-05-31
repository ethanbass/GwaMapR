# GwaMapR 0.2.0

* Added `plot_gwas` function to plot gene arrow diagrams for multiple genes (after clumping GWAS results). The old `plot_gwas` function has been renamed as `plot_gwas_single`.
* Added `calculate_maf` function to calculating the minor allele frequency of a specified SNP from the BED file.
* Added `extract_genes` function to export gene or protein sequences for selected genes as FASTA files.
* Added `get_genes` function to return list of genes associated with a specified SNP.
* Added `get_peaks` function to return a list of genes associated with a phenotype based on results from `gemma` and a `GFF3` annotation file.
* Added `read_gff` function to read a GFF3 file into R.
* Added `get_snps` function to extract raw SNP data from a bed file for specified SNPs.
* Fixed `title` argument in `plot_manhattan`.
