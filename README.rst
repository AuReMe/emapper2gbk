.. image:: https://img.shields.io/pypi/v/emapper2gbk.svg
	:target: https://pypi.org/project/emapper2gbk

.. image:: https://img.shields.io/github/license/AuReMe/emapper_to_gbk.svg
	:target: https://github.com/AuReMe/emapper_to_gbk/blob/master/LICENSE

.. image:: https://github.com/AuReMe/emapper_to_gbk/workflows/Python%20package/badge.svg
    :target: https://github.com/AuReMe/emapper_to_gbk/actions

.. image:: https://img.shields.io/badge/doi-10.7554/eLife.61968-blueviolet.svg
	:target: https://doi.org/10.7554/eLife.61968

emapper2gbk: creation of genbank files from Eggnog-mapper annotation outputs
============================================================================

Starting from fasta and `Eggnog-mapper <http://eggnog-mapper.embl.de/>`__ annotation files, build a gbk file that is suitable for metabolic network reconstruction with `Pathway Tools <http://bioinformatics.ai.sri.com/ptools/>`__. Adds the GO terms and EC numbers annotations in the genbank file.
There are two main modes:

* **genes mode**: suitable when a list of isolated genes/proteins have been annotated with Eggnog-mapper, typically the gene catalogue of a metagenome.

* **genomes mode**: usually when focusing on a single organism, with a ``.gff`` file. The creation of genbanks can be performed in parallel by providing directories (with matching names for genomes, proteomes and annotation files) as inputs.

**If you use emapper2gbk, please cite**

Belcour* A, Frioux* C, Aite M, Bretaudeau A, Hildebrand F, Siegel A. Metage2Metabo, microbiota-scale metabolic complementarity for the identification of key species. eLife 2020;9:e61968 `https://doi.org/10.7554/eLife.61968 <https://doi.org/10.7554/eLife.61968>`_ .

.. contents:: Table of contents
   :backlinks: top
   :local:

Main inputs
-----------

emapper2gbk genes
~~~~~~~~~~~~~~~~~

For each annotated list of genes, inputs are:

* a nucleotide fasta file containing the CDS sequence of each genes or a folder containing multiple nucleotide fasta files.
* the translated sequences in amino-acids in fasta or a folder containing the corresponding protein sequences to the nucleotide sequences (must be the same name).
* the annotation file obtained after Eggnog-mapper annotation (usually ``xxx.emapper.annotation``) or a folder with multiple annotation files (must be the same name as nucleotide fasta file and ends with '.tsv' extension).

In addition, as optional files:

* the name of the considered organism (can be "bacteria" or "metagenome"), a full lineage of the organism (such as `Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia,Escherichia coli`) or a file with organisms names (matching the genomes names).
* the merge option to merge genes into fake contigs.
* the number of available cores for multiprocessing (when working on multiple genomes).
* a go-basic file of GO ontology (if not given, emapper2gbk will download a copy and use it).

Example:

Input with files:

.. code-block:: text

    nucleotide_sequences.fna
    protein_sequence.faa
    annotation.emapper.annotation

Input with folders:

.. code-block:: text

    nucleotide_sequences
    ├── gene_list_1.fna
    ├── gene_list_2.fna
    protein_sequence
    ├── gene_list_1.faa
    ├── gene_list_2.faa
    annotation
    ├── gene_list_1.tsv
    ├── gene_list_2.tsv

.. image:: pictures/emapper2gbk_genes.svg

To work the ID of the genes in the nucleic fasta file (``-fn``) must be the same than the ID of the proteins in the protein fasta file (``-fp``) and in the annotation file (``-a``).

emapper2gbk genomes
~~~~~~~~~~~~~~~~~~~

For each genomes, inputs are:

* a nucleotide fasta file containing the sequence of each contigs/chromosomes for the genome or a folder containing multiple nucleotide fasta files.
* the proteome corresponding to the genome or a folder containing the corresponding protein sequences to the nucleotide sequences (having the same name as the nucleotides files).
* the GFF file corresponding to the genome or a folder containing multiple GFF files (each GFF files must have the same name as the corresponding nucleotide files).
* the annotation file obtained after Eggnog-mapper annotation (usually ``xxx.emapper.annotation``) or a folder with multiple annotation files (must be the same name as nucleotide fasta file and ends with '.tsv' extension)

In addition, as optional files:

* the name of the considered organism (can be "bacteria"), a full lineage of the organism (such as ``Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia,Escherichia coli``) or a file with organisms names (matching the genomes names).
* the number of available cores for multiprocessing (when working on multiple genomes).
* a go-basic file of GO ontology (if not given, emapper2gbk will download a copy and use it).

Example:

Input with files:

.. code-block:: text

    nucleotide_sequences.fna
    protein_sequence.faa
    annotation.emapper.annotation
    genome.gff

Input with folders:

.. code-block:: text

    nucleotide_sequences
    ├── genome_1.fna
    ├── genome_2.fna
    protein_sequence
    ├── genome_1.faa
    ├── genome_2.faa
    annotation
    ├── genome_1.tsv
    ├── genome_2.tsv
    gff
    ├── genome_1.gff
    ├── genome_2.gff

.. image:: pictures/emapper2gbk_genomes.svg

The ID in the chromosome/contigs/scaffolds fasta file (``-fn``) must correspond to region in the gff file (``-g``).
Then the genes in the region will be found and the child CDS associated to the genes wil be extracted.
The CDS ID must be the same than the ID in the protein fasta file (``-fp``) and the ID in the eggnog-mapper annotation file (``-a``).


By default emapper2gbk searches for inheritance between genes and CDS in the GFF files.
A gene feature is required and the CDS feature must have the gene feature as a parent, like in this example:

.. code-block:: text

    ##gff file
    region_1	RefSeq	region	1	12642	.	+	.	ID=region_1
    region_1	RefSeq	gene	1	2445	.	-	.	ID=gene_1
    region_1	RefSeq	CDS	1	2445	.	-	0	ID=cds_1;Parent=gene_1

Depending on which field (CDS, mRNA or gene) of the gff is associated with the proteome IDs in faa file, the gff-type (``-gt``) option can take into account these 3 parameters (``CDS``, ``mRNA``, ``gene``).
The tool also takes into account particular gff formats (Gmove and eggnog) and the gff-type option (``-gt``) can take these 2 parameters (``gmove``, ``eggnog``).

**CDS gff type**

For example, some GFF files can be formatted differently with only CDS (such as in `Prodigal <https://github.com/hyattpd/Prodigal>`__ or `Prokka <https://github.com/tseemann/prokka>`__ GFF), it is possible to use them with ``-gt CDS`` (case sensitive).
Here is an example of the format accepted by this command (with ID cds_1 being the same as the one in the faa and eggnogg-mapper files):

.. code-block:: text

    ##gff file
    region_1	RefSeq	CDS	1	2445	.	-	0	ID=cds_1

**mRNA gff type**

The ``-gt mRNA`` option (case sensitive) is to be used in case the protein identifiers in the faa file match the identifiers in the "mRNA" field of the gff.
Here is an example of the format accepted by this command (with ID cds_1 being the same as the one in the faa and eggnogg-mapper files):

.. code-block:: text

    ##gff file
    region_1	RefSeq	mRNA	1	2445	.	-	0	ID=cds_1


It is useful for gff formats containing multiple "CDS" fields associated with 1 gene and/or 1 mRNA and a single sequence in the faa file.

**gene gff type**

The ``-gt gene`` option is to be used in case the protein identifiers in the faa file match the identifiers in the "gene" field of the gff.
Here is an example of the format accepted by this command (with ID cds_1 being the same as the one in the faa and eggnogg-mapper files):

.. code-block:: text

    ##gff file
    region_1	RefSeq	gene	1	2445	.	-	0	ID=cds_1


**Gmove gff type**

The tool handle GFF from `Gmove <https://www.genoscope.cns.fr/gmove/>`__ (with ``-gt gmove``) with the following format:

.. code-block:: text

    ##gff file
    region_1	Gmove	mRNA	1	2445	.	+	.	ID=mRNA_gene_1;Name=mRNA_gene_1
    region_1	Gmove	CDS	1	2445	.	-	0	Parent=mRNA_gene_1

For gmove, the proteins in the faa and eggnogg-mapper files will be prefixed with ``prot_`` (like ``prot_gene_1`` for ``mRNA_gene_1``). Emapper2gbk should be able to handle these differences.

**EggNog gff type**

It is also possible to use the GFF created by eggnog-mapper (if a fasta genome was given as input to eggnog-mapper) with ``-gt eggnog``.
An example of such use can be seen in the `test folder <https://github.com/AuReMe/emapper2gbk/tree/master/tests/test_data/data_eggnog>`__.

Taxonomic information
~~~~~~~~~~~~~~~~~~~~~

There is 3 possible ways to give taxonomic information to emapper2gbk:

* ``-n "Scientific name"``: using only the -n option, it is possible to give a scientific name of an organism (compliant with the NCBI Taxonomy database). This name will be queried against the EBI to extract taxonomic information.

* ``-n "Kingdom;Order;Class;Family;Genus;Species`` --ete`: adding the ``--ete`` parameter will change how `-n` works, it will then expect a full lineage (compliant with NCBI Taxonomy database, such as ``Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia,Escherichia coli``). This will be parsed by the ete3 package to extract the taxonomic information.

* ``-nf taxonomic_information.tsv``: for multiple genomes, it is possible to use the option ``-nf``. It expects a tsv file with a first column containing name of the input files and a second column with the scientific name (or lineage) of the associated organism. An example (`organism_names.tsv <https://github.com/AuReMe/emapper2gbk/blob/main/tests/test_data/organism_names.tsv>`__) is present in the test folder.

Dependencies and installation
-----------------------------

Dependencies
~~~~~~~~~~~~

All are described in ``requirements.txt`` and can be installed with ``pip install -r requirements.txt``.

* biopython
* ete
* gffutils
* pandas
* pronto
* requests

Install
~~~~~~~

* From this cloned repository

.. code-block:: sh

    pip install -r requirements.txt
    pip install .

* From Pypi

.. code-block:: sh

    pip install emapper2gbk

Usage
-----

Convert GFF, fastas, annotation table and species name into Genbank.

.. code-block:: sh

    usage: emapper2gbk [-h] [-v] {genes,genomes} ...

    Starting from fasta and Eggnog-mapper annotation files, build a gbk file that is suitable for metabolic network reconstruction with Pathway Tools. Adds the GO terms and EC numbers annotations in the genbank file.

    Two modes:
    - genomes (one genome/proteome/gff/annot file --> one gbk).
    - genes with the annotation of the full gene catalogue and fasta files (nucleic and protein) corresponding to list of genes.

    Examples:

    * Genomic - single mode

    emapper2gbk genomes -fn genome.fna -fp proteome.faa -gff genome.gff -n "Escherichia coli" -o coli.gbk -a eggnog_annotation.tsv [-go go-basic.obo]

    * Genomic - multiple mode, "bacteria" as default name

    emapper2gbk genes -fn genome_dir/ -fp proteome_dir/ -n metagenome -o gbk_dir/ -a eggnog_annotation_dir/ [-go go-basic.obo]

    * Genomic - multiple mode, tsv file for organism names

    emapper2gbk genes -fn genome_dir/ -fp proteome_dir/ -nf matching_genome_orgnames.tsv -o gbk_dir/ -a eggnog_annotation_dir/ [-go go-basic.obo]

    * Metagenomic

    emapper2gbk genes -fn genome_dir/ -fp proteome_dir/ -o gbk_dir/ -a gene_cat_ggnog_annotation.tsv --one-annot-file [-go go-basic.obo]

    You can give the GO ontology as an input to the program, it will be otherwise downloaded during the run. You can download it here: http://purl.obolibrary.org/obo/go/go-basic.obo .
    The program requests the NCBI database to retrieve taxonomic information of the organism. However, if the organism is "bacteria" or "metagenome", the taxonomic information will not have to be retrieved online.
    Hence, if you need to run the program from a cluster with no internet access, it is possible for a "bacteria" or "metagenome" organism, and by providing the GO-basic.obo file.
    For specific help on each subcommand use: emapper2gbk {cmd} --help

    optional arguments:
    -h, --help       show this help message and exit
    -v, --version    show program's version number and exit

    subcommands:
    valid subcommands:

    {genes,genomes}
        genes          genes mode : 1-n annot, 1-n faa, 1-n fna (gene sequences) --> 1 gbk
        genomes        genomes mode: 1-n contig/chromosome fasta, 1-n protein fasta, 1-n GFF, 1-n annot --> 1 gbk


* Genomes mode

  * Usage

    .. code-block:: sh

        usage: emapper2gbk genomes [-h] -fn FASTANUCLEIC -fp FASTAPROT -o OUPUT_DIR -g GFF [-gt GFF_TYPE] [-nf NAMEFILE] [-n NAME] -a ANNOTATION [-c CPU] [-go GOBASIC] [-q] [--keep-gff-annotation] [--ete]

        Build a gbk file for each genome with an annotation file for each

        options:
        -h, --help            show this help message and exit
        -fn FASTANUCLEIC, --fastanucleic FASTANUCLEIC
                                fna file or directory
        -fp FASTAPROT, --fastaprot FASTAPROT
                                faa file or directory
        -o OUPUT_DIR, --out OUPUT_DIR
                                output directory/file path
        -g GFF, --gff GFF     gff file or directory
        -gt GFF_TYPE, --gff-type GFF_TYPE
                                gff type, by default emapper2gbk search for CDS with gene as Parent in the GFF, but by using the '-gt cds_only' option emapper2gbk will only use the CDS information from the genome, by using '-gt gmove' emapper2gbk will use mRNA to find CDS, by
                                using 'eggnog' emapper2gbk will use the output ifles of eggnog-mapper
        -nf NAMEFILE, --namefile NAMEFILE
                                organism/genome name (col 2) associated to genome file basenames (col 1). Default = 'metagenome' for metagenomic and 'cellular organisms' for genomic
        -n NAME, --name NAME  organism/genome name in quotes
        -a ANNOTATION, --annotation ANNOTATION
                                eggnog annotation file or directory
        -c CPU, --cpu CPU     cpu number for metagenomic mode or genome mode using input directories
        -go GOBASIC, --gobasic GOBASIC
                                go ontology, GOBASIC is either the name of an existing file containing the GO Ontology or the name of the file that will be created by emapper2gbk containing the GO Ontology
        -q, --quiet           quiet mode, only warning, errors logged into console
        --keep-gff-annotation
                                Copy the annotation from the GFF (product) into the genbank output file.
        --ete                 Use ete3 NCBITaxa instead of query on the EBI Taxonomy Database for taxonomic ID assignation (useful if there is no internet access, except that ete3 NCBITaxa database must have been downloaded before).

  * Examples

    * Genomic - single mode

    .. code:: sh

      emapper2gbk genomes -fn genome.fna -fp proteome.faa -gff genome.gff -n "Escherichia coli" -o coli.gbk -a eggnog_annotation.tsv [-go go-basic.obo]

    * Genomic - multiple mode, "bacteria" as default name

* genes mode

  * Usage

    .. code-block:: sh

        usage: emapper2gbk genes [-h] -fn FASTANUCLEIC -fp FASTAPROT -o OUPUT_DIR -a ANNOTATION [-c CPU] [-n NAME] [-nf NAMEFILE] [-go GOBASIC] [--merge MERGE] [-q] [--ete]

        Use the annotation of a complete gene catalogue and build gbk files for each set of genes (fna) and proteins (faa) from input directories

        options:
        -h, --help            show this help message and exit
        -fn FASTANUCLEIC, --fastanucleic FASTANUCLEIC
                                fna file or directory
        -fp FASTAPROT, --fastaprot FASTAPROT
                                faa file or directory
        -o OUPUT_DIR, --out OUPUT_DIR
                                output directory/file path
        -a ANNOTATION, --annotation ANNOTATION
                                eggnog annotation file or directory
        -c CPU, --cpu CPU     cpu number for metagenomic mode or genome mode using input directories
        -n NAME, --name NAME  organism/genome name in quotes
        -nf NAMEFILE, --namefile NAMEFILE
                                organism/genome name (col 2) associated to genome file basenames (col 1). Default = 'metagenome' for metagenomic and 'cellular organisms' for genomic
        -go GOBASIC, --gobasic GOBASIC
                                go ontology, GOBASIC is either the name of an existing file containing the GO Ontology or the name of the file that will be created by emapper2gbk containing the GO Ontology
        --merge MERGE         Number of gene sequences to merge into fake contig from a same file in the genbank file.
        -q, --quiet           quiet mode, only warning, errors logged into console
        --ete                 Use ete3 NCBITaxa instead of query on the EBI Taxonomy Database for taxonomic ID assignation (useful if there is no internet access, except that ete3 NCBITaxa database must have been downloaded before).

  * Example

    .. code:: sh

      emapper2gbk genes -fn genome_dir/ -fp proteome_dir/ -o gbk_dir/ -a gene_cat_ggnog_annotation.tsv [-go go-basic.obo]

License
-------

This project is licensed under the GNU LGPL-3.0-or-later - see the `LICENSE file <https://github.com/AuReMe/emapper2gbk/blob/main/LICENSE>`__ for details.