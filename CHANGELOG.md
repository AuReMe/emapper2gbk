# Changelog

# emapper2gbk v0.1.0 (2021-06-29)

This release modifies heavily `emapper2gbk`. The `genomic` and `metagenomic` subcommands have been removed. From this release, emapper2gbk is divided in two subcommands: `genes` and `genomes`.
The `genes` subcommand can be used for a catalogue  of genes (with a fasta file for gene nucleic acid sequence, a fasta file for the protein sequence and the eggnog-mapper annotaiton file).
The `genomes` subcommand can be user for genome (with a fasta file for genome sequences, one for protein sequences, the eggnog-mapper anntoation file and a GFF file).
The two subcommands can be used on single file or on folder containing multiple files.

The scripts files of emapper2gbk have been renamed according to these changes.

## Add

* New subcommand `emapper2gbk genes`: subcommand to use for catalogue of genes (with three inputs: fasta of nucleic sequences of genes, fasta of amino acids sequences of proteins and the eggnog-mapper annotation file).
* New subcommand `emapper2gbk genomes`: subcommand to use for genomes (with four inputs: fasta of nucleic sequences of genomes (chromosomes sequences), fasta of amino acids sequences of proteins, the eggnog-mapper annotation file and a GFF file).
* New option `-gt` for `emapper2gbk genomes`: `-gt cds_only` if the GFF file contains only CDS information.
* New option `--ete` for `emapper2gbk genomes`: when searching for taxonomic ID instead of requesting the EBI, use the package `ete3`.
*  New option `--keep-gff-annotation`  for `emapper2gbk genomes`: keep the `product` field of the CDS in the GFF file into the genbank file.
* Add dbxref (Kegg, Bigg, Pfam and CAZ) to genbank.

## Modify

* Replace `-fg` argument of fasta nucleic sequences (for genes or genomes) by `-fn`.
* Replace GPL license by LGPL license.
* Refactoring of the code of emapper2gbk.
* Update readme (add pictures, badges and modify text about new subcommands).

## Fix

* Numerous issues with the format of the output of new version of eggnog-mapper.
* Issue with Preferred_name.
* Issue with GO Terms.

## Remove

* Remove  `genomic` and `metagenomic` subcommands.
* Remove unused files.
  
# emapper_to_gbk v0.0.7 (2020-09-08)

## Fix

* Issue with replacement of Bio.Alphabet in Biopython 1.78. This version is compatible with both versions: old versions (with Bio.Alphabet) and new versions (without Bio.Alphabet).

# emapper_to_gbk v0.0.6 (2020-06-24)

## Fix

* Issue with version of setup.cfg with PyPI.

# emapper_to_gbk v0.0.5 (2020-06-24)

## Fix

* Issue in long_description of setup.py with PyPI.

# emapper_to_gbk v0.0.4 (2020-06-24)

## Fix

* Issue in long_description of setup.py with PyPI.

# emapper_to_gbk v0.0.3 (2020-06-24)

## Fix

* Issue in long_description of setup.py with PyPI.

# emapper_to_gbk v0.0.2 (2020-06-24)

## Fix

* Issue in long_description of setup.py with PyPI.

# emapper_to_gbk v0.0.1 (2020-06-24)

**First release of emapper2gbk**

## Features

* `emapper2gbk` entrypoint
* `emapper2gbk genomic` to create genbank from eggnog-mapper results and single organism (with or without GFF)
* `emapper2gbk metagenomic` to create genbank from eggnog-mapper results and gene catalogue of a metagenome (with or without GFF)

## Integration

* Use Github Actions for CI
* Release on PyPI