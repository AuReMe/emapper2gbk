# Copyright (C) 2019-2025 Clémence Frioux & Arnaud Belcour - Inria Dyliss - Pleiade
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

"""Console script for emapper2gbk."""
import argparse
import logging
import os
import sys
import time

from argparse import RawTextHelpFormatter

from emapper2gbk import __version__ as VERSION
from emapper2gbk.emapper2gbk import gbk_creation
from emapper2gbk.utils import is_valid_file, is_valid_path

LICENSE = """Copyright (C) Pleiade and Dyliss Inria projects\n
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.\n

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.\n

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>\n
"""

MESSAGE = """
Starting from fasta and Eggnog-mapper annotation files, build a gbk file that is suitable for metabolic network reconstruction with Pathway Tools. Adds the GO terms and EC numbers annotations in the genbank file.\n
Two modes:
- genomes (one genome/proteome/gff/annot file --> one gbk).
- genes with the annotation of the full gene catalogue and fasta files (nucleic and protein) corresponding to list of genes. \n

Examples: \n
* One genome of "Escherichia coli" \n
emapper2gbk genomes -fn genome.fna -fp proteome.faa -gff genome.gff -n "Escherichia coli" -o coli.gbk -a eggnog_annotation.tsv [-go go-basic.obo] \n
* Multiple genomes \n
emapper2gbk genes -fn genome_dir/ -fp proteome_dir/ -n metagenome -o gbk_dir/ -a eggnog_annotation_dir/ [-go go-basic.obo] \n
* One genes list \n
emapper2gbk genes -fn genes.fna -fp genes.faa -o genes.gbk -a genes.emapper.annotation [-go go-basic.obo] \n
* Multiple genes list \n
emapper2gbk genes -fn genes_dir/ -fp proteomes_dir/ -nf matching_genome_orgnames.tsv -o gbk_dir/ -a eggnog_annotation_dir/ [-go go-basic.obo] \n
* Multiple genes list with one annotation file \n
emapper2gbk genes -fn genes_dir/ -fp proteomes_dir/ -o gbk_dir/ -a gene_cat_ggnog_annotation.tsv [-go go-basic.obo]
\n

You can give the GO ontology as an input to the program, it will be otherwise downloaded during the run. You can download it here: http://purl.obolibrary.org/obo/go/go-basic.obo .
The program requests the NCBI database to retrieve taxonomic information of the organism. However, if the organism is "bacteria", "metagenome", "archaea" or "eukaryota", the taxonomic information will not have to be retrieved online.
Hence, if you need to run the program from a cluster with no internet access, it is possible for a "bacteria", "metagenome", "archaea" or "eukaryota" organism, and by providing the GO-basic.obo file.
"""

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def cli():
    """Console script for emapper2gbk."""
    start_time = time.time()
    parser = argparse.ArgumentParser(
        "emapper2gbk",
        description=MESSAGE + " For specific help on each subcommand use: emapper2gbk {cmd} --help", formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s " + VERSION + "\n" + LICENSE)

    # parent parsers
    parent_parser_q = argparse.ArgumentParser(add_help=False)
    parent_parser_q.add_argument(
        "-q",
        "--quiet",
        dest="quiet",
        help="quiet mode, only warning, errors logged into console",
        required=False,
        action="store_true",
        default=None,
    )
    parent_parser_c = argparse.ArgumentParser(add_help=False)
    parent_parser_c.add_argument(
        "-c",
        "--cpu",
        help="cpu number for metagenomic mode or genome mode using input directories",
        required=False,
        type=int,
        default=1
    )
    parent_parser_o = argparse.ArgumentParser(add_help=False)
    parent_parser_o.add_argument(
        "-o",
        "--out",
        dest="out",
        required=True,
        help="output directory/file path",
        metavar="OUPUT_DIR"
    )
    parent_parser_faa = argparse.ArgumentParser(add_help=False)
    parent_parser_faa.add_argument(
        "-fp",
        "--fastaprot",
        help="faa file or directory",
        required=True,
        type=str
    )
    parent_parser_fna = argparse.ArgumentParser(add_help=False)
    parent_parser_fna.add_argument(
        "-fn",
        "--fastanucleic",
        help="fna file or directory",
        required=True,
        type=str
    )
    parent_parser_gff = argparse.ArgumentParser(add_help=False)
    parent_parser_gff.add_argument(
        "-g",
        "--gff",
        help="gff file or directory",
        required=True,
        type=str
    )
    parent_parser_gff_type = argparse.ArgumentParser(add_help=False)
    parent_parser_gff_type.add_argument(
        "-gt",
        "--gff-type",
        help="gff type, by default emapper2gbk search for CDS with gene as Parent in the GFF, but by using the '-gt cds_only' option emapper2gbk will only use the CDS information from the genome, by using '-gt gmove' emapper2gbk will use mRNA to find CDS, by using 'eggnog' emapper2gbk will use the output ifles of eggnog-mapper",
        required=False,
        type=str
    )
    parent_parser_ann = argparse.ArgumentParser(add_help=False)
    parent_parser_ann.add_argument(
        "-a",
        "--annotation",
        help="eggnog annotation file or directory",
        required=True,
        type=str
    )
    parent_parser_go = argparse.ArgumentParser(add_help=False)
    parent_parser_go.add_argument(
        "-go",
        "--gobasic",
        help="go ontology, GOBASIC is either the name of an existing file containing the GO Ontology or the name of the file that will be created by emapper2gbk containing the GO Ontology",
        required=False,
        default=None,
        type=str
    )
    parent_parser_name = argparse.ArgumentParser(add_help=False)
    parent_parser_name.add_argument(
        "-n",
        "--name",
        help="organism/genome name in quotes",
        required=False,
        # default="Bacteria",
        type=str
    )
    parent_parser_namef = argparse.ArgumentParser(add_help=False)
    parent_parser_namef.add_argument(
        "-nf",
        "--namefile",
        help="organism/genome name (col 2) associated to genome file basenames (col 1). Default = 'metagenome' for metagenomic and 'cellular organisms' for genomic",
        required=False,
        type=str)
    parent_parser_merge = argparse.ArgumentParser(add_help=False)
    parent_parser_merge.add_argument(
        "--merge",
        dest="merge",
        help="Number of gene sequences to merge into fake contig from a same file in the genbank file.",
        required=False,
        type=int,
        default=None
    )
    parent_parser_keep_gff_annot = argparse.ArgumentParser(add_help=False)
    parent_parser_keep_gff_annot.add_argument(
        "--keep-gff-annotation",
        dest="keep_gff_annotation",
        help="Copy the annotation from the GFF (product) into the genbank output file.",
        required=False,
        action="store_true",
        default=None,
    )
    parent_parser_ete = argparse.ArgumentParser(add_help=False)
    parent_parser_ete.add_argument(
        "--ete",
        dest="ete",
        help="Use ete4 NCBITaxa instead of query on the EBI Taxonomy Database for taxonomic ID assignation (useful if there is no internet access, except that ete4 NCBITaxa database must have been downloaded before).",
        required=False,
        action="store_true",
        default=None,
    )

   # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest="cmd")
    genes_parser = subparsers.add_parser(
        "genes",
        help="genes mode : 1-n annot, 1-n faa, 1-n fna (gene sequences) --> 1 gbk",
        parents=[
            parent_parser_fna, parent_parser_faa, parent_parser_o,
            parent_parser_ann, parent_parser_c, parent_parser_name, parent_parser_namef,
            parent_parser_go, parent_parser_merge, parent_parser_q, parent_parser_ete
        ],
        description=
        "Use the annotation of a complete gene catalogue and build gbk files for each set of genes (fna) and proteins (faa) from input directories",
        allow_abbrev=False
    )
    genomes_parser = subparsers.add_parser(
        "genomes",
        help="genomes mode: 1-n contig/chromosome fasta, 1-n protein fasta, 1-n GFF, 1-n annot --> 1 gbk",
        parents=[
            parent_parser_fna, parent_parser_faa, parent_parser_o, parent_parser_gff, parent_parser_gff_type,
            parent_parser_namef, parent_parser_name, parent_parser_ann,
            parent_parser_c, parent_parser_go, parent_parser_q, parent_parser_keep_gff_annot,
            parent_parser_ete
        ],
        description=
        "Build a gbk file for each genome with an annotation file for each",
        allow_abbrev=False
    )

    args = parser.parse_args()

    # If no argument print the help.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # test writing in out_directory if a subcommand is given else print version and help
    if args.cmd:
        if not is_valid_path(args.out):
            logger.critical("Impossible to access/create output directory/file")
            sys.exit(1)
    else:
        logger.info("emapper2gbk " + VERSION + "\n" + LICENSE)
        parser.print_help()
        sys.exit()

    # # add logger in file
    formatter = logging.Formatter('%(message)s')

    # set up the default console logger
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)

    if args.quiet:
        console_handler.setLevel(logging.WARNING)
    logger.addHandler(console_handler)

    # check go-basic file
    if args.gobasic and not is_valid_file(args.gobasic):
        logger.critical(f"No Go-basic file available, it will be download by emapper2gbk.")

    args = parser.parse_args()

    # check the given names
    if args.namefile and args.name:
        logger.warning("You should either use a --name or --namefile, not both. Will consider the file only.")
        orgnames = args.namefile
    elif args.namefile:
        orgnames = args.namefile
    elif args.name:
        orgnames = args.name
    else:
        if args.cmd == "genomes":
            orgnames = "cellular organisms"
            logger.warning("The default organism name 'cellular organisms' is used.")
        if args.cmd == "genes":
            orgnames = "metagenome"
            logger.warning("The default organism name 'metagenome' is used.")

    # Check name.
    if args.namefile:
        if os.path.isfile(args.fastanucleic) and os.path.isfile(args.fastaprot):
            logger.error("Tabulated file for organisms name should not be used for single runs of genomic mode. Will use the --name argument or the default 'metagenome'for metagenomic or 'cellular organisms' for genomics name if None")
            if args.name:
                orgnames = args.name
            else:
                if args.cmd == "genomes":
                    orgnames = "cellular organisms"
                    logger.warning("The default organism name 'cellular organisms' is used.")
                if args.cmd == "genes":
                    orgnames = "metagenome"
                    logger.warning("The default organism name 'metagenome' is used.")

    if args.cmd == "genomes":
        if not args.gff_type:
            gff_type = 'default'
        else:
            gff_type = args.gff_type
        gbk_creation(nucleic_fasta=args.fastanucleic, protein_fasta=args.fastaprot, annot=args.annotation, gff=args.gff, gff_type=gff_type,
                        org=orgnames, output_path=args.out, gobasic=args.gobasic, cpu=args.cpu, keep_gff_annot=args.keep_gff_annotation,
                        ete_option=args.ete)

    elif args.cmd == "genes":
        gbk_creation(nucleic_fasta=args.fastanucleic, protein_fasta=args.fastaprot, annot=args.annotation, org=orgnames,
                        output_path=args.out, gobasic=args.gobasic, cpu=args.cpu, merge_genes_fake_contig=args.merge,
                        ete_option=args.ete)


    logger.info("--- Total runtime %.2f seconds ---" % (time.time() - start_time))

if __name__ == "__main__":
    sys.exit(cli())  
