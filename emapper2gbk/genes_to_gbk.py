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

import logging

from Bio import SeqFeature as sf
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict
from typing import Union

from emapper2gbk.utils import create_cds_feature, check_valid_path, create_GO_namespaces_alternatives, read_annotation, \
                                create_taxonomic_data, create_taxonomic_data_ete, get_basename, record_info

logger = logging.getLogger(__name__)

"""
Description:
Using fasta files (scaffold/chromosme/contig file, protein file), annotation tsv file from eggnog and the species name
this script writes a genbank file with EC number and Go annotations.
The species name needs to be compatible with the taxonomy of the EBI.
Informations need a good formating:
gene ID should be correctly written (like XXX_001 and no XXX_1 if you got more thant 100 genes).
Currently when there is multiple GO terms/EC the script split them when they are separated by ";" or by "," like GO:0006979;GO:0020037;GO:0004601,
if you use another separator add to the re.split(',|;').
Other informations can be added by adding a dictionary with gene ID as key and the information
as value and adapt the condition used for the others annotations (EC, Go term).
"""


def faa_to_gbk(nucleic_fasta:str, protein_fasta:str, annot:Union[str, dict],
                org:str, output_path:str, gobasic:Union[None, str, dict],
                merge_genes_fake_contig:int, ete_option:bool):
    """ Create genbank file from nucleic and protein fasta plus eggnog mapper annotation file.

    Args:
        nucleic_fasta (str): nucleic fasta file
        protein_fasta (str): protein fasta file
        annot (str): annotation file or dictionary
        org (str): organims name or mapping file
        output_path (str): output file or directory
        gobasic (str): path to go-basic.obo file or dictionary
        merge_genes_fake_contig (int): merge genes into fake contig. The int associted to merge is the number of genes per fake contigs.
        ete_option (bool): to use ete4 NCBITaxa database for taxonomic ID assignation instead of request on the EBI taxonomy database.
    """
    check_valid_path([nucleic_fasta, protein_fasta])

    genome_id = get_basename(nucleic_fasta)

    # Dictionary with gene id as key and nucleic sequence as value.
    gene_nucleic_seqs = OrderedDict()

    for record in SeqIO.parse(nucleic_fasta, "fasta"):
        gene_nucleic_seqs[record.id] = record.seq

    # Dictionary with gene id as key and protein sequence as value.
    gene_protein_seqs = OrderedDict()

    for record in SeqIO.parse(protein_fasta, "fasta"):
        protein_id = record.id
        if protein_id.isnumeric():
            protein_id = f"gene_{protein_id}"
        gene_protein_seqs[protein_id] = record.seq

    # Create a taxonomy dictionary querying the EBI.
    if ete_option:
        species_informations = create_taxonomic_data_ete(org)
    else:
        species_informations = create_taxonomic_data(org)
    if species_informations is None:
        return False

    # Read the eggnog tsv file containing GO terms and EC associated with gene name.
    # if metagenomic mode, annotation is already read and given as a dict
    if not type(annot) is dict:
        annot = dict(read_annotation(annot))

    # Query Gene Ontology to extract namespaces and alternative IDs.
    # go_namespaces: Dictionary GO id as term and GO namespace as value.
    # go_alternatives: Dictionary GO id as term and GO alternatives id as value.
    if gobasic:
        if not type(gobasic[0]) is dict and not type(gobasic[1]) is dict:
            go_namespaces, go_alternatives = create_GO_namespaces_alternatives(gobasic)
        else:
            go_namespaces, go_alternatives = gobasic
    else:
        go_namespaces, go_alternatives = create_GO_namespaces_alternatives()

    logger.info('Assembling Genbank informations for ' + genome_id)

    # Create fake contig by merging genes.
    if merge_genes_fake_contig:
        create_genbank_fake_contig(gene_nucleic_seqs, gene_protein_seqs, annot, go_namespaces, go_alternatives, output_path, species_informations, merge_genes_fake_contig)
    else:
        create_genbank(gene_nucleic_seqs, gene_protein_seqs, annot, go_namespaces, go_alternatives, output_path, species_informations)

    return True

def create_genbank(gene_nucleic_seqs, gene_protein_seqs, annot,
                    go_namespaces, go_alternatives, output_path,
                    species_informations):
    """ Create genbank file from nucleic and protein fasta plus eggnog mapper annotation file.

    Args:
        gene_nucleic_seqs (dict): dictionary of nucleic sequences (key: sequence id, value: sequence)
        gene_protein_seqs (dict): dictionary of protein sequences (key: sequence id, value: sequence)
        annot (dict): dictionary of eggnog-ammper annotation (key: gene_id, value: ['GOs','EC', 'Preferred_name'])
        go_namespaces (dict): dictionary of GO terms namespace (key: GO Term ID, value: namespace associated to GO Term)
        go_alternatives (dict): dictionary of GO terms alternatives ID (key: GO Term ID, value: alternatives GO Term associated to GO Term)
        output_path (str): output file or directory
        species_informations (dict): dictionary containing information about species
    """
    # All SeqRecord objects will be stored in a list and then give to the SeqIO writer to create the genbank.
    records = []

    # Iterate through each contig/gene.
    for gene_nucleic_id in sorted(gene_nucleic_seqs):
        # Create a SeqRecord object using gene information.
        record = record_info(gene_nucleic_id, gene_nucleic_seqs[gene_nucleic_id], species_informations)

        # If id is numeric, change it
        if gene_nucleic_id.isnumeric():
            id_gene = f"gene_{gene_nucleic_id}"
        elif "|" in gene_nucleic_id:
            id_gene = gene_nucleic_id.split("|")[1]
        else:
            id_gene = gene_nucleic_id
        start_position = 1
        end_position = len(gene_nucleic_seqs[gene_nucleic_id])
        strand = 0
        new_feature_gene = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                            end_position,
                                                            strand),
                                                            type="gene")
        new_feature_gene.qualifiers['locus_tag'] = id_gene

        # Add gene information to contig record.
        record.features.append(new_feature_gene)

        new_cds_feature = create_cds_feature(id_gene, start_position, end_position, strand, annot, go_namespaces, go_alternatives, gene_protein_seqs)
        new_cds_feature.qualifiers['locus_tag'] = id_gene

        # Add CDS information to contig record
        record.features.append(new_cds_feature)

        records.append(record)

    # Create Genbank with the list of SeqRecord.
    SeqIO.write(records, output_path, 'genbank')


def create_genbank_fake_contig(gene_nucleic_seqs, gene_protein_seqs, annot,
                            go_namespaces, go_alternatives, output_path,
                            species_informations, merge_genes_fake_contig):
    """ Create genbank file from nucleic and protein fasta plus eggnog mapper annotation file.

    Args:
        gene_nucleic_seqs (dict): dictionary of nucleic sequences (key: sequence id, value: sequence)
        gene_protein_seqs (dict): dictionary of protein sequences (key: sequence id, value: sequence)
        annot (dict): dictionary of eggnog-ammper annotation (key: gene_id, value: ['GOs','EC', 'Preferred_name'])
        go_namespaces (dict): dictionary of GO terms namespace (key: GO Term ID, value: namespace associated to GO Term)
        go_alternatives (dict): dictionary of GO terms alternatives ID (key: GO Term ID, value: alternatives GO Term associated to GO Term)
        output_path (str): output file or directory
        species_informations (dict): dictionary containing information about species
        merge_genes_fake_contig (int): merge genes into fake contig. The int associted to merge is the number of genes per fake contigs.
    """
    fake_contigs = {}
    # Iterate through each contig/gene.
    maximal_contig = int(len(list(gene_nucleic_seqs.keys())) / merge_genes_fake_contig)

    maximal_contig_str_stize = len(str(maximal_contig))

    def homogonize_id(contig_id, max_id_length):
        old_contig_number= contig_id.split('_')[-1]
        if len(old_contig_number) < max_id_length:
            nb_diff = max_id_length - len(old_contig_number)
            new_contig_id = '_'.join(contig_id.split('_')[:-1]) + '_' + str(nb_diff*[0]) + str(old_contig_number)
        else:
            new_contig_id = contig_id

        return new_contig_id

    fake_contigs.update({homogonize_id('fake_contig_'+str(index), maximal_contig_str_stize):list(gene_nucleic_seqs.keys())[x:x+merge_genes_fake_contig] for index, x in enumerate(range(0, len(gene_nucleic_seqs), merge_genes_fake_contig))})

    # All SeqRecord objects will be stored in a list and then give to the SeqIO writer to create the genbank.
    records = []
    for contig in fake_contigs:
        contig_seq = Seq('').join([gene_nucleic_seqs[gene] for gene in fake_contigs[contig]])

        # Create a SeqRecord object using gene information.
        record = record_info(contig, contig_seq, species_informations)

        gene_position = 0
        for id_gene in fake_contigs[contig]:
            start_position = gene_position
            gene_position += len(gene_nucleic_seqs[id_gene])
            end_position = gene_position
            strand = 0
            new_feature_gene = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                                end_position,
                                                                strand),
                                                                type="gene")
            new_feature_gene.qualifiers['locus_tag'] = id_gene

            # Add gene information to contig record.
            record.features.append(new_feature_gene)

            new_cds_feature = create_cds_feature(id_gene, start_position, end_position, strand, annot, go_namespaces, go_alternatives, gene_protein_seqs)
            new_cds_feature.qualifiers['locus_tag'] = id_gene

            # Add CDS information to contig record
            record.features.append(new_cds_feature)

        records.append(record)
    SeqIO.write(records, output_path, 'genbank')
