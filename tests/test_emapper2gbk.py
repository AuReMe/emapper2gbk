#!/usr/bin/env python

"""Tests for emapper2gbk package."""

import os
import shutil
import subprocess

from Bio import SeqIO
from emapper2gbk.emapper2gbk import gbk_creation


FAA_INPUT = os.path.join('test_data', 'betaox_genes.faa')
FNA_INPUT = os.path.join('test_data', 'betaox_genes.fna')
NUMERIC_FNA_INPUT = os.path.join('test_data', 'numeric', 'betaox_genes.fna')
NUMERIC_FAA_INPUT = os.path.join('test_data', 'numeric', 'betaox_genes.faa')

GENOME_FNA_DIR = os.path.join('test_data', 'fna_genomes_mode')
GENOME_FNA_INPUT = os.path.join('test_data', 'betaox_genomes.fna')
GENOME_FAA_INPUT = os.path.join('test_data', 'betaox_genomes.faa')
GENOME_GFF_INPUT = os.path.join('test_data', 'betaox_genomes.gff')
NUMERIC_ANNOT_GFF = os.path.join('test_data', 'numeric', 'betaox_genomes.gff')

ANNOT_INPUT = os.path.join('test_data', 'betaox_annotation.tsv')
NUMERIC_ANNOT_INPUT = os.path.join('test_data', 'numeric', 'betaox_annotation.tsv')
GENOME_ANNOT_INPUT = os.path.join('test_data', 'betaox_annotation_genomes.tsv')
ANNOT_INPUT_V2 = os.path.join('test_data', 'betaox_v2.emapper.annotations')
ANNOT_DIR = os.path.join('test_data', 'ann')
ANNOT_DIR_GENOME = os.path.join('test_data', 'ann_genomes')

FAA_DIR = os.path.join('test_data', 'faa')
FAA_DIR_GENOMES = os.path.join('test_data', 'faa_genomes')
FNA_DIR = os.path.join('test_data', 'fna')
GFF_DIR = os.path.join('test_data', 'gff')

ORG_NAME = 'Escherichia coli'
ORG_NAME_BACT = 'bacteria'
ORG_NAME_ARCH = 'archaea'
ORG_NAME_EUK = 'eukaryota'
ORG_NAME_META = 'metagenome'
ORG_NAME_CELL = 'cellular organisms'
ORG_FULL_TAX = 'cellular organisms;Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli'
ORG_FILE = os.path.join('test_data', 'organism_names.tsv')
GO_FILE = 'go-basic.obo'

EXPECTED_GBK_NO_GFF = os.path.join('test_data', 'betaox_no_gff.gbk')
EXPECTED_GBK_WITH_GFF = os.path.join('test_data', 'betaox_from_gff.gbk')
EXPECTED_GBK_NO_GFF_MERGED = os.path.join('test_data', 'betaox_no_gff_merged.gbk')
EXPECTED_GBK_NUMERIC = os.path.join('test_data', 'numeric', 'betaox.gbk')
EXPECTED_GBK_GFF_NUMERIC = os.path.join('test_data', 'numeric', 'betaox_genomes.gbk')

GMOVE_ANNOT = os.path.join('test_data', 'data_gmove', 'betaox_v2.emapper.annotations')
GMOVE_FNA = os.path.join('test_data', 'data_gmove', 'betaox_genomes.fna')
GMOVE_FAA = os.path.join('test_data', 'data_gmove', 'betaox_genomes.faa')
GMOVE_GFF = os.path.join('test_data', 'data_gmove', 'betaox_genomes.gff')

EGGNOG_ANNOT = os.path.join('test_data', 'data_eggnog', 'out.emapper.annotations')
EGGNOG_FNA = os.path.join('test_data', 'data_eggnog', 'queries.fasta')
EGGNOG_FAA = os.path.join('test_data', 'data_eggnog', 'out.emapper.genepred.fasta')
EGGNOG_GFF = os.path.join('test_data', 'data_eggnog', 'out.emapper.genepred.gff')
EGGNOG_GBK = os.path.join('test_data', 'data_eggnog', 'eggnog_expected.gbk')

TYPE_CDS_GFF = os.path.join('test_data', 'data_gff_type_accession', 'betaox_genomes_CDS.gff')
TYPE_MRNA_GFF = os.path.join('test_data', 'data_gff_type_accession', 'betaox_genomes_mRNA.gff')
TYPE_GENE_GFF = os.path.join('test_data', 'data_gff_type_accession', 'betaox_genomes_gene.gff')

FASTA_ONLY_GENOME_FNA_DIR = os.path.join('test_data', 'fasta_only', 'fna_genomes_mode')
FASTA_ONLY_GENOME_FAA_DIR = os.path.join('test_data', 'fasta_only', 'faa_genomes')
FASTA_ONLY_GENOME_ANN_DIR = os.path.join('test_data', 'fasta_only', 'ann_genomes')

ANNOTATIONS_TYPES = ['go_function', 'go_process', 'go_component', 'EC_number', 'locus', 'db_xref']

ANNOTATIONS_BY_GENOME = {'gene1781':{'go_component':['GO:0005575', 'GO:0005623', 'GO:0005886',
                                                    'GO:0016020', 'GO:0044464', 'GO:0071944'],
                                    'go_function':['GO:0003674', 'GO:0003824', 'GO:0015645',
                                                    'GO:0016405', 'GO:0016874', 'GO:0016877',
                                                    'GO:0016878', 'GO:0031956'],
                                    'go_process':[],
                                    'EC_number':[]
                                    },
                        'gene1887':{'go_component':['GO:0005575', 'GO:0005622', 'GO:0005623',
                                                    'GO:0005737', 'GO:0005829', 'GO:0005886',
                                                    'GO:0009898', 'GO:0016020', 'GO:0044424',
                                                    'GO:0044425', 'GO:0044444', 'GO:0044459',
                                                    'GO:0044464', 'GO:0071944', 'GO:0098552',
                                                    'GO:0098562'],
                                    'go_function':['GO:0003674', 'GO:0003824', 'GO:0004467',
                                                    'GO:0005488', 'GO:0005504', 'GO:0008289',
                                                    'GO:0015645', 'GO:0016405', 'GO:0016874',
                                                    'GO:0016877', 'GO:0016878', 'GO:0031406',
                                                    'GO:0033293', 'GO:0036041', 'GO:0036094',
                                                    'GO:0043167', 'GO:0043168', 'GO:0043177',
                                                    'GO:0070538'],
                                    'go_process':['GO:0001676', 'GO:0006082', 'GO:0006139',
                                                'GO:0006163', 'GO:0006629', 'GO:0006631',
                                                'GO:0006635', 'GO:0006637', 'GO:0006644',
                                                'GO:0006725', 'GO:0006732', 'GO:0006753',
                                                'GO:0006790', 'GO:0006793', 'GO:0006796',
                                                'GO:0006807', 'GO:0008150', 'GO:0008152',
                                                'GO:0008610', 'GO:0008654', 'GO:0009056',
                                                'GO:0009058', 'GO:0009062', 'GO:0009117',
                                                'GO:0009150', 'GO:0009259', 'GO:0009314',
                                                'GO:0009411', 'GO:0009416', 'GO:0009628',
                                                'GO:0009987', 'GO:0016042', 'GO:0016054',
                                                'GO:0019395', 'GO:0019637', 'GO:0019693',
                                                'GO:0019752', 'GO:0030258', 'GO:0032787',
                                                'GO:0033865', 'GO:0033875', 'GO:0034032',
                                                'GO:0034440', 'GO:0034641', 'GO:0035383',
                                                'GO:0043436', 'GO:0043603', 'GO:0044237',
                                                'GO:0044238', 'GO:0044242', 'GO:0044248',
                                                'GO:0044249', 'GO:0044255', 'GO:0044281',
                                                'GO:0044282', 'GO:0046395', 'GO:0046483',
                                                'GO:0050896', 'GO:0051186', 'GO:0055086',
                                                'GO:0055114', 'GO:0071704', 'GO:0072329',
                                                'GO:0072521', 'GO:0090407', 'GO:1901135',
                                                'GO:1901360', 'GO:1901564', 'GO:1901575',
                                                'GO:1901576'],
                                    'EC_number':['6.2.1.3']
                                    },
                        'gene2441':{'go_function':['GO:0003674', 'GO:0003824', 'GO:0003857',
                                                    'GO:0004165', 'GO:0004300', 'GO:0008691',
                                                    'GO:0008692', 'GO:0016491', 'GO:0016614',
                                                    'GO:0016616', 'GO:0016829', 'GO:0016835',
                                                    'GO:0016836', 'GO:0016853', 'GO:0016854',
                                                    'GO:0016856', 'GO:0016860', 'GO:0016863'],
                                    'go_process':['GO:0006082', 'GO:0006629', 'GO:0006631',
                                                    'GO:0006635', 'GO:0006725', 'GO:0006805',
                                                    'GO:0008150', 'GO:0008152', 'GO:0009056',
                                                    'GO:0009062', 'GO:0009404', 'GO:0009407',
                                                    'GO:0009410', 'GO:0009636', 'GO:0009850',
                                                    'GO:0009852', 'GO:0009987', 'GO:0010124',
                                                    'GO:0010817', 'GO:0016042', 'GO:0016054',
                                                    'GO:0019395', 'GO:0019439', 'GO:0019748',
                                                    'GO:0019752', 'GO:0030258', 'GO:0032787',
                                                    'GO:0034440', 'GO:0042178', 'GO:0042221',
                                                    'GO:0042445', 'GO:0042447', 'GO:0042537',
                                                    'GO:0043436', 'GO:0044237', 'GO:0044238',
                                                    'GO:0044242', 'GO:0044248', 'GO:0044255',
                                                    'GO:0044281', 'GO:0044282', 'GO:0046395',
                                                    'GO:0050896', 'GO:0051716', 'GO:0055114',
                                                    'GO:0065007', 'GO:0065008', 'GO:0070887',
                                                    'GO:0071466', 'GO:0071704', 'GO:0072329',
                                                    'GO:0098754', 'GO:1901360', 'GO:1901361',
                                                    'GO:1901575'],
                                    'EC_number':['1.1.1.157',
                                                '1.1.1.35',
                                                '4.2.1.17',
                                                '5.1.2.3',
                                                '5.3.3.8'
                                                ]
                                    },
                        'gene3987':{'go_component':['GO:0005575', 'GO:0005622', 'GO:0005623',
                                                    'GO:0005737', 'GO:0044424', 'GO:0044464'],
                                    'go_function':['GO:0003674', 'GO:0003824', 'GO:0003988',
                                                    'GO:0016408', 'GO:0016740', 'GO:0016746',
                                                    'GO:0016747'],
                                    'go_process':['GO:0006082', 'GO:0006629', 'GO:0006631',
                                                'GO:0006635', 'GO:0008150', 'GO:0008152',
                                                'GO:0009056', 'GO:0009062', 'GO:0009987',
                                                'GO:0016042', 'GO:0016054', 'GO:0019395',
                                                'GO:0019752', 'GO:0030258', 'GO:0032787',
                                                'GO:0034440', 'GO:0043436', 'GO:0044237',
                                                'GO:0044238', 'GO:0044242', 'GO:0044248',
                                                'GO:0044255', 'GO:0044281', 'GO:0044282',
                                                'GO:0046395', 'GO:0055114', 'GO:0071704',
                                                'GO:0072329', 'GO:1901575'],
                                    'EC_number':['2.3.1.16']
                                    },
                        'gene3988':{'go_component':[],
                                    'go_function':['GO:0003674', 'GO:0003824', 'GO:0003857',
                                                'GO:0004165', 'GO:0004300', 'GO:0008691',
                                                'GO:0008692', 'GO:0016491', 'GO:0016614',
                                                'GO:0016616', 'GO:0016829', 'GO:0016835',
                                                'GO:0016836', 'GO:0016853', 'GO:0016854',
                                                'GO:0016856', 'GO:0016860', 'GO:0016863'],
                                    'go_process':['GO:0006082', 'GO:0006629', 'GO:0006631',
                                                'GO:0006635', 'GO:0006725', 'GO:0006805',
                                                'GO:0008150', 'GO:0008152', 'GO:0009056',
                                                'GO:0009062', 'GO:0009404', 'GO:0009407',
                                                'GO:0009410', 'GO:0009636', 'GO:0009850',
                                                'GO:0009852', 'GO:0009987', 'GO:0010124',
                                                'GO:0010817', 'GO:0016042', 'GO:0016054',
                                                'GO:0019395', 'GO:0019439', 'GO:0019748',
                                                'GO:0019752', 'GO:0030258', 'GO:0032787',
                                                'GO:0034440', 'GO:0042178', 'GO:0042221',
                                                'GO:0042445', 'GO:0042447', 'GO:0042537',
                                                'GO:0043436', 'GO:0044237', 'GO:0044238',
                                                'GO:0044242', 'GO:0044248', 'GO:0044255',
                                                'GO:0044281', 'GO:0044282', 'GO:0046395',
                                                'GO:0050896', 'GO:0051716', 'GO:0055114',
                                                'GO:0065007', 'GO:0065008', 'GO:0070887',
                                                'GO:0071466', 'GO:0071704', 'GO:0072329',
                                                'GO:0098754', 'GO:1901360', 'GO:1901361',
                                                'GO:1901575'],
                                    'EC_number':['1.1.1.157',
                                                '1.1.1.35',
                                                '4.2.1.17',
                                                '5.1.2.3',
                                                '5.3.3.8'
                                                ]
                                    }
                        }


def compare_two_gbks(expected_gbk:str, tested_gbk:str):
    """Compare the annotations of 2 genbank files.

    Args:
        expected_gbk (str): path to expected gbk
        tested_gbk (str): path to the second gbk
    """
    loaded_gbk = SeqIO.to_dict(SeqIO.parse(expected_gbk, "genbank"))
    loaded_test = SeqIO.to_dict(SeqIO.parse(tested_gbk, "genbank"))

    assert set(loaded_gbk.keys()) == set(loaded_test.keys())

    annotations_expected = {i:None
                            for i in loaded_gbk.keys()}

    annotations_test = {i:None
                        for i in loaded_test.keys()}

    for gene in annotations_expected:
        for index in range(0, len(loaded_gbk[gene].features)):
            if loaded_gbk[gene].features[index].type=='CDS':
                annotations_expected[gene] = loaded_gbk[gene].features[index].qualifiers

    for gene in annotations_test:
        for index in range(0, len(loaded_test[gene].features)):
            if loaded_test[gene].features[index].type=='CDS':
                annotations_test[gene] = loaded_test[gene].features[index].qualifiers

    for gene in annotations_expected:
        for qualifier in annotations_expected[gene]:
            if qualifier in ANNOTATIONS_TYPES:
                assert set(annotations_expected[gene][qualifier]) == set(annotations_test[gene][qualifier])

    return


def check_gbks_from_dir_genes_mode(gbk_dir):
    """Check if annotations in each gbk file are consistent with the expected ones.

    Args:
        gbk_dir (str): path to gbk directory
    """
    for gbk in os.listdir(gbk_dir):
        gbk_path = os.path.join(gbk_dir, gbk)
        loaded_gbk = SeqIO.to_dict(SeqIO.parse(gbk_path, "genbank"))
        annotations = {i:None
                        for i in loaded_gbk.keys()}
        for gene in annotations:
            for index in range(0, len(loaded_gbk[gene].features)):
                if loaded_gbk[gene].features[index].type=='CDS':
                    annotations[gene] = loaded_gbk[gene].features[index].qualifiers
                    # check annotations
                    for ann in ANNOTATIONS_TYPES:
                        if ann in annotations[gene]:
                            assert set(annotations[gene][ann]) == set(ANNOTATIONS_BY_GENOME[gene][ann])

    return


def check_gbks_from_dir_genome_mode(gbk_dir):
    """Check if annotations in each gbk file are consistent with the expected ones.

    Args:
        gbk_dir (str): path to gbk directory
    """
    expected_genomes = {}
    for expected_gene in ANNOTATIONS_BY_GENOME:
        expected_genomes['cds-'+expected_gene] = ANNOTATIONS_BY_GENOME[expected_gene]
    for gbk in os.listdir(gbk_dir):
        gbk_path = os.path.join(gbk_dir, gbk)
        for record in SeqIO.parse(gbk_path, "genbank"):
            for feature in record.features:
                if feature.type == 'CDS':
                    gene = feature.qualifiers['locus_tag'][0]
                    qualifier_annotations = {qualifier: feature.qualifiers[qualifier] for qualifier in feature.qualifiers}
                    annotations = {}
                    annotations[gene] = qualifier_annotations
                    for ann in ANNOTATIONS_TYPES:
                        if ann in annotations[gene]:
                            assert set(annotations[gene][ann]) == set(expected_genomes[gene][ann])

    return


def check_gbks_from_dir_genome_mode_fasta_only(gbk_dir):
    """Check if annotations in each gbk file are consistent with the expected ones.

    Args:
        gbk_dir (str): path to gbk directory
    """
    expected_genomes = {}
    for expected_gene in ANNOTATIONS_BY_GENOME:
        expected_genomes['NC_000913_3_'+expected_gene] = ANNOTATIONS_BY_GENOME[expected_gene]
    for gbk in os.listdir(gbk_dir):
        gbk_path = os.path.join(gbk_dir, gbk)
        for record in SeqIO.parse(gbk_path, "genbank"):
            for feature in record.features:
                if feature.type == 'CDS':
                    gene = feature.qualifiers['locus_tag'][0]
                    qualifier_annotations = {qualifier: feature.qualifiers[qualifier] for qualifier in feature.qualifiers}
                    annotations = {}
                    annotations[gene] = qualifier_annotations
                    for ann in ANNOTATIONS_TYPES:
                        if ann in annotations[gene]:
                            assert set(annotations[gene][ann]) == set(expected_genomes[gene][ann])

    return


def test_gbk_genes_mode_test():
    """Test genes mode with file as input.
    """
    print("*** Test genes mode with file as input ***")
    gbk_test = 'test_no_gff.gbk'
    gbk_creation(nucleic_fasta=FNA_INPUT,
                protein_fasta=FAA_INPUT,
                annot=ANNOT_INPUT,
                org=ORG_NAME,
                output_path=gbk_test,
                gff=None,
                gobasic=GO_FILE)

    compare_two_gbks(EXPECTED_GBK_NO_GFF, gbk_test)
    os.remove(gbk_test)

    return

def test_gbk_genes_mode_lineage_test():
    """Test genes mode with file as input and full lineage taxonomy.
    """
    print("*** Test genes mode with file as input and full lineage taxonomy ***")
    gbk_test = 'test_no_gff.gbk'
    gbk_creation(nucleic_fasta=FNA_INPUT,
                protein_fasta=FAA_INPUT,
                annot=ANNOT_INPUT,
                org=ORG_FULL_TAX,
                output_path=gbk_test,
                gff=None,
                gobasic=GO_FILE)

    compare_two_gbks(EXPECTED_GBK_NO_GFF, gbk_test)
    os.remove(gbk_test)

    return


def test_gbk_gene_mode_test_cli():
    """Test genes mode  with file as input.
    """
    gbk_test = 'test_no_gff.gbk'

    print("*** Test genes mode with file as input with cli***")
    subprocess.call(['emapper2gbk', 'genes', '-fn', FNA_INPUT, '-fp', FAA_INPUT,
                        '-a', ANNOT_INPUT, '-o', gbk_test, '-go', GO_FILE, '-n', ORG_NAME_BACT])

    compare_two_gbks(EXPECTED_GBK_NO_GFF, gbk_test)
    os.remove(gbk_test)

    return


def test_gbk_genes_mode_annot_v2():
    """Test genes mode with file as input with v2 eggnog mapper annotation file.
    """
    print("*** Test genes mode with file as input with v2 eggnog mapper annotation file***")
    gbk_test = 'test_no_gff.gbk'
    gbk_creation(nucleic_fasta=FNA_INPUT,
                protein_fasta=FAA_INPUT,
                annot=ANNOT_INPUT_V2,
                org=ORG_NAME_ARCH,
                output_path=gbk_test,
                gff=None,
                gobasic=GO_FILE)

    compare_two_gbks(EXPECTED_GBK_NO_GFF, gbk_test)
    os.remove(gbk_test)

    return


def test_gbk_genes_mode_annot_v2_cli():
    """Test genes mode with file as input with v2 eggnog mapper annotation file with cli.
    """
    gbk_test = 'test_no_gff.gbk'
    print("*** Test genes mode with file as input with v2 eggnog mapper annotation file with cli***")
    subprocess.call(['emapper2gbk', 'genes', '-fn', FNA_INPUT, '-fp', FAA_INPUT,
                        '-a', ANNOT_INPUT_V2, '-o', gbk_test, '-go', GO_FILE, '-n', ORG_NAME_BACT])

    compare_two_gbks(EXPECTED_GBK_NO_GFF, gbk_test)
    os.remove(gbk_test)

    return


def test_gbk_genomes_mode_test():
    """Test genomes mode with file as input.
    """
    print("*** Test genomes mode with file as input ***")
    gbk_test = 'test_gff.gbk'

    gbk_creation(nucleic_fasta=GENOME_FNA_INPUT,
                protein_fasta=GENOME_FAA_INPUT,
                annot=GENOME_ANNOT_INPUT,
                gff=GENOME_GFF_INPUT,
                org=ORG_NAME_ARCH,
                output_path=gbk_test,
                gobasic=GO_FILE)

    compare_two_gbks(EXPECTED_GBK_WITH_GFF, gbk_test)
    os.remove(gbk_test)

    return

def test_gbk_genomes_mode_lineage_test():
    """Test genomes mode with file as input and full lineage taxonomy.
    """
    print("*** Test genomes mode with file as input and full lineage taxonomy ***")
    gbk_test = 'test_gff.gbk'

    gbk_creation(nucleic_fasta=GENOME_FNA_INPUT,
                protein_fasta=GENOME_FAA_INPUT,
                annot=GENOME_ANNOT_INPUT,
                gff=GENOME_GFF_INPUT,
                org=ORG_FULL_TAX,
                output_path=gbk_test,
                gobasic=GO_FILE)

    compare_two_gbks(EXPECTED_GBK_WITH_GFF, gbk_test)
    os.remove(gbk_test)

    return


def test_gbk_genomes_mode_test_cli():
    """Test genomes mode with file as input with cli
    """
    gbk_test = 'test_gff.gbk'
    print("*** Test genomes mode with file as input with cli***")
    subprocess.call(['emapper2gbk', 'genomes', '-fn', GENOME_FNA_INPUT, '-fp', GENOME_FAA_INPUT,
                        '-a', GENOME_ANNOT_INPUT, '-g', GENOME_GFF_INPUT, '-o', gbk_test, '-go', GO_FILE,
                        '-n', ORG_NAME_EUK])

    compare_two_gbks(EXPECTED_GBK_WITH_GFF, gbk_test)
    os.remove(gbk_test)

    return


def test_gbk_genomes_mode_CDS_test():
    """Test genomes mode with file as input.
    """
    print("*** Test genomes mode with file as input and using CDS option ***")
    gbk_test = 'test_gff.gbk'

    gbk_creation(nucleic_fasta=GENOME_FNA_INPUT,
                protein_fasta=GENOME_FAA_INPUT,
                annot=GENOME_ANNOT_INPUT,
                gff=TYPE_CDS_GFF,
                gff_type='CDS',
                org=ORG_NAME_EUK,
                output_path=gbk_test,
                gobasic=GO_FILE)

    compare_two_gbks(EXPECTED_GBK_WITH_GFF, gbk_test)
    os.remove(gbk_test)

    return


def test_gbk_genomes_mode_mRNA_test():
    """Test genomes mode with file as input.
    """
    print("*** Test genomes mode with file as input and using mRNA option ***")
    gbk_test = 'test_gff.gbk'

    gbk_creation(nucleic_fasta=GENOME_FNA_INPUT,
                protein_fasta=GENOME_FAA_INPUT,
                annot=GENOME_ANNOT_INPUT,
                gff=TYPE_MRNA_GFF,
                gff_type='mRNA',
                org=ORG_NAME_EUK,
                output_path=gbk_test,
                gobasic=GO_FILE)

    compare_two_gbks(EXPECTED_GBK_WITH_GFF, gbk_test)
    os.remove(gbk_test)

    return


def test_gbk_genomes_mode_gene_test():
    """Test genomes mode with file as input.
    """
    print("*** Test genomes mode with file as input and using gene option ***")
    gbk_test = 'test_gff.gbk'

    gbk_creation(nucleic_fasta=GENOME_FNA_INPUT,
                protein_fasta=GENOME_FAA_INPUT,
                annot=GENOME_ANNOT_INPUT,
                gff=TYPE_GENE_GFF,
                gff_type='gene',
                org=ORG_NAME_EUK,
                output_path=gbk_test,
                gobasic=GO_FILE)

    compare_two_gbks(EXPECTED_GBK_WITH_GFF, gbk_test)
    os.remove(gbk_test)

    return

def test_gbk_genomes_mode_keep_gff_annot_test():
    """Test genomes mode with file as input.
    """
    print("*** Test genomes mode with file as input and using keep_gff_annotation option***")
    gbk_test = 'test_gff.gbk'

    gbk_creation(nucleic_fasta=GENOME_FNA_INPUT,
                protein_fasta=GENOME_FAA_INPUT,
                annot=GENOME_ANNOT_INPUT,
                gff=GENOME_GFF_INPUT,
                org=ORG_NAME_META,
                output_path=gbk_test,
                gobasic=GO_FILE,
                keep_gff_annot=True)

    compare_two_gbks(EXPECTED_GBK_WITH_GFF, gbk_test)
    os.remove(gbk_test)

    return


def test_gbk_genomes_mode_folder():
    """Test genomes mode with folders as input
    """
    print("*** Test genomes mode with folders as input ***")
    gbk_dir_test = 'gbk_mg'
    gbk_creation(nucleic_fasta=GENOME_FNA_DIR,
                protein_fasta=FAA_DIR_GENOMES,
                annot=ANNOT_DIR_GENOME,
                org=ORG_FILE,
                gff=GFF_DIR,
                output_path=gbk_dir_test,
                gobasic=GO_FILE)

    check_gbks_from_dir_genome_mode(gbk_dir_test)
    shutil.rmtree(gbk_dir_test)

    return


def test_gbk_genomes_mode_folder_cli():
    """Test genomes mode with folders as input with cli
    """
    gbk_dir_test = 'gbk_mg'

    print("*** Test genomes mode with folders as input with cli ***")
    subprocess.call(['emapper2gbk', 'genomes', '-fn', GENOME_FNA_DIR, '-fp', FAA_DIR_GENOMES,
                        '-a', ANNOT_DIR_GENOME, '-g', GFF_DIR, '-o', gbk_dir_test, '-go', GO_FILE,
                        '-nf', ORG_FILE])

    check_gbks_from_dir_genome_mode(gbk_dir_test)
    shutil.rmtree(gbk_dir_test)

    return


def test_gbk_genes_mode_folder():
    """Test genes mode with folders as input
    """
    print("*** Test genes mode with folders as input ***")
    gbk_dir_test = 'gbk_g'

    gbk_creation(nucleic_fasta=FNA_DIR,
                protein_fasta=FAA_DIR,
                annot=ANNOT_DIR,
                org=ORG_FILE,
                gff=None,
                output_path=gbk_dir_test,
                gobasic=GO_FILE)

    check_gbks_from_dir_genes_mode(gbk_dir_test)
    shutil.rmtree(gbk_dir_test)

    return


def test_gbk_genes_mode_folder_cli():
    """Test genes mode with folders as input with cli
    """
    gbk_dir_test = 'gbk_g'

    print("*** Test genes mode with folders as input with cli ***")
    subprocess.call(['emapper2gbk', 'genes', '-fn', FNA_DIR, '-fp', FAA_DIR,
                        '-a', ANNOT_DIR, '-o', gbk_dir_test, '-go', GO_FILE,
                        '-nf', ORG_FILE])

    check_gbks_from_dir_genes_mode(gbk_dir_test)
    shutil.rmtree(gbk_dir_test)

    return


def test_gbk_genes_mode_folder_one_annot_file():
    """Test genes mode with folders as input with one annotation file
    """
    print("*** Test genes mode with folders as input with one annotation file ***")
    gbk_dir_test = 'gbk_mg'

    gbk_creation(nucleic_fasta=FNA_DIR,
                protein_fasta=FAA_DIR,
                annot=ANNOT_INPUT,
                org=ORG_FILE,
                gff=None,
                output_path=gbk_dir_test,
                gobasic=GO_FILE)

    check_gbks_from_dir_genes_mode(gbk_dir_test)
    shutil.rmtree(gbk_dir_test)

    return


def test_gbk_genes_mode_folder_one_annot_file_cli():
    """Test genes mode with folders as input with one annotation file with cli
    """
    gbk_dir_test = 'gbk_mg'

    print("*** Test genes mode with folders as input with one annotation file with cli ***")
    subprocess.call(['emapper2gbk', 'genes', '-fn', FNA_DIR, '-fp', FAA_DIR,
                        '-a', ANNOT_INPUT, '-o', gbk_dir_test, '-go', GO_FILE,
                        '-nf', ORG_FILE])

    check_gbks_from_dir_genes_mode(gbk_dir_test)
    shutil.rmtree(gbk_dir_test)

    return


def test_gbk_genes_mode_merge_fake_contig():
    """Test genes mode with files as input and merge genes in fake contig
    """
    print("*** Test genes mode with files as input and merge genes in fake contig ***")
    gbk_test = 'test_merged.gbk'
    gbk_creation(nucleic_fasta=FNA_INPUT,
                protein_fasta=FAA_INPUT,
                annot=ANNOT_INPUT,
                org=ORG_NAME_BACT,
                output_path=gbk_test,
                gff=None,
                gobasic=GO_FILE,
                merge_genes_fake_contig=5)

    # Dictionary with gene id as key and nucleic sequence as value.
    expected_gene_nucleic_seqs = dict()
    for record in SeqIO.parse(FNA_INPUT, 'fasta'):
        expected_gene_nucleic_seqs[record.id] = record.seq

    # Fasta sequence form genbank
    found_gene_nucleic_seq = dict()
    for record in SeqIO.parse(gbk_test, 'genbank'):
        for feature in record.features:
            if feature.type == 'gene':
                expected_gene_nucleic_seqs[feature.qualifiers['locus_tag'][0]] = feature.extract(record.seq)

    # Check if the sequences from the genbank are corrects.
    for gene_id in expected_gene_nucleic_seqs:
        assert expected_gene_nucleic_seqs[gene_id] == expected_gene_nucleic_seqs[gene_id]

    compare_two_gbks(EXPECTED_GBK_NO_GFF_MERGED, gbk_test)
    os.remove(gbk_test)

    return


def test_gbk_genomes_mode_gmove_test():
    """Test genomes mode with file as input.
    """
    print("*** Test genomes mode with file as input ***")
    gbk_test = 'test_gff.gbk'

    gbk_creation(nucleic_fasta=GMOVE_FNA,
                protein_fasta=GMOVE_FAA,
                annot=GMOVE_ANNOT,
                gff=GMOVE_GFF,
                gff_type='gmove',
                org=ORG_NAME_BACT,
                output_path=gbk_test,
                gobasic=GO_FILE)

    compare_two_gbks(EXPECTED_GBK_WITH_GFF, gbk_test)
    os.remove(gbk_test)

    return


def test_gbk_genomes_mode_gmove_test_cli():
    """Test genomes mode with file as input with cli
    """
    gbk_test = 'test_gff.gbk'
    print("*** Test genomes mode with file as input with cli***")
    subprocess.call(['emapper2gbk', 'genomes', '-fn', GMOVE_FNA, '-fp', GMOVE_FAA,
                        '-a', GMOVE_ANNOT, '-g', GMOVE_GFF, '-o', gbk_test, '-go', GO_FILE,
                        '-n', ORG_NAME_BACT, '-gt', 'gmove'])

    compare_two_gbks(EXPECTED_GBK_WITH_GFF, gbk_test)
    os.remove(gbk_test)

    return

def test_gbk_genomes_mode_eggnog_test():
    """Test genomes mode with file as input.
    """
    print("*** Test genomes mode eggnog with file as input ***")
    gbk_test = 'test_gff.gbk'

    gbk_creation(nucleic_fasta=EGGNOG_FNA,
                protein_fasta=EGGNOG_FAA,
                annot=EGGNOG_ANNOT,
                gff=EGGNOG_GFF,
                gff_type='eggnog',
                org=ORG_NAME_BACT,
                output_path=gbk_test,
                gobasic=GO_FILE)

    compare_two_gbks(EGGNOG_GBK, gbk_test)
    os.remove(gbk_test)

    return


def test_gbk_genomes_mode_eggnog_test_cli():
    """Test genomes mode with file as input with cli
    """
    gbk_test = 'test_gff.gbk'
    print("*** Test genomes mode eggnog with file as input with cli***")
    subprocess.call(['emapper2gbk', 'genomes', '-fn', EGGNOG_FNA, '-fp', EGGNOG_FAA,
                        '-a', EGGNOG_ANNOT, '-g', EGGNOG_GFF, '-o', gbk_test, '-go', GO_FILE,
                        '-n', ORG_NAME_BACT, '-gt', 'eggnog'])

    compare_two_gbks(EGGNOG_GBK, gbk_test)
    os.remove(gbk_test)

    return

def test_gbk_genes_mode_numeric_test():
    """Test genes mode with file as input.
    """
    print("*** Test genes mode with file as input ***")
    gbk_test = 'test_no_gff.gbk'

    gbk_creation(nucleic_fasta=NUMERIC_FNA_INPUT,
                protein_fasta=NUMERIC_FAA_INPUT,
                annot=NUMERIC_ANNOT_INPUT,
                org=ORG_NAME,
                output_path=gbk_test,
                gff=None,
                gobasic=GO_FILE)

    compare_two_gbks(EXPECTED_GBK_NUMERIC, gbk_test)
    os.remove(gbk_test)

    return

def test_gbk_genes_mode_numeric_ete_test():
    """Test genes mode with file as input.
    """
    print("*** Test genes mode with file as input ***")
    gbk_test = 'test_no_gff.gbk'

    gbk_creation(nucleic_fasta=NUMERIC_FNA_INPUT,
                protein_fasta=NUMERIC_FAA_INPUT,
                annot=NUMERIC_ANNOT_INPUT,
                org=ORG_NAME,
                output_path=gbk_test,
                gff=None,
                gobasic=GO_FILE)

    compare_two_gbks(EXPECTED_GBK_NUMERIC, gbk_test)
    os.remove(gbk_test)

    return

def test_gbk_genomes_numeric_test():
    """Test genomes mode with file as input.
    """
    print("*** Test genomes mode eggnog with file as input ***")
    gbk_test = 'test_gff.gbk'

    gbk_creation(nucleic_fasta=GENOME_FNA_INPUT,
                protein_fasta=NUMERIC_FAA_INPUT,
                annot=NUMERIC_ANNOT_INPUT,
                gff=NUMERIC_ANNOT_GFF,
                org=ORG_NAME_BACT,
                output_path=gbk_test,
                ete_option=True,
                gobasic=GO_FILE)

    compare_two_gbks(EXPECTED_GBK_GFF_NUMERIC, gbk_test)
    os.remove(gbk_test)

    return


def test_gbk_genomes_mode_folder_fasta_only():
    """Test genomes mode with folders as input
    """
    print("*** Test genomes mode with folders as input ***")
    gbk_dir_test = 'gbk_mg'
    gbk_creation(nucleic_fasta=FASTA_ONLY_GENOME_FNA_DIR,
                protein_fasta=FASTA_ONLY_GENOME_FAA_DIR,
                annot=FASTA_ONLY_GENOME_ANN_DIR,
                org=ORG_FILE,
                gff='fasta-only',
                gff_type='fasta-only',
                output_path=gbk_dir_test,
                gobasic=GO_FILE)

    check_gbks_from_dir_genome_mode_fasta_only(gbk_dir_test)
    shutil.rmtree(gbk_dir_test)

    return


def test_gbk_genomes_mode_folder_fasta_only_cli():
    """Test genomes mode with folders as input
    """
    print("*** Test genomes mode with folders as input ***")
    gbk_dir_test = 'gbk_mg'
    subprocess.call(['emapper2gbk', 'genomes', '-fn', FASTA_ONLY_GENOME_FNA_DIR, '-fp', FASTA_ONLY_GENOME_FAA_DIR,
                        '-a', FASTA_ONLY_GENOME_ANN_DIR, '-g', 'fasta-only', '-gt', 'fasta-only', '-o', gbk_dir_test, '-go', GO_FILE,
                        '-n', ORG_FILE])

    check_gbks_from_dir_genome_mode_fasta_only(gbk_dir_test)
    shutil.rmtree(gbk_dir_test)

    return

if __name__ == "__main__":
    test_gbk_genes_mode_test()
    test_gbk_genes_mode_lineage_test()
    test_gbk_gene_mode_test_cli()
    test_gbk_genes_mode_annot_v2()
    test_gbk_genes_mode_annot_v2_cli()
    test_gbk_genomes_mode_test()
    test_gbk_genomes_mode_lineage_test()
    test_gbk_genomes_mode_test_cli()
    test_gbk_genomes_mode_folder()
    test_gbk_genomes_mode_folder_cli()
    test_gbk_genes_mode_folder()
    test_gbk_genes_mode_folder_cli()
    test_gbk_genes_mode_folder_one_annot_file()
    test_gbk_genes_mode_folder_one_annot_file_cli()
    test_gbk_genes_mode_merge_fake_contig()
    test_gbk_genomes_mode_gmove_test()
    test_gbk_genomes_mode_gmove_test_cli()
    test_gbk_genes_mode_numeric_test()
    test_gbk_genomes_numeric_test()
    test_gbk_genomes_mode_CDS_test()
    test_gbk_genomes_mode_mRNA_test()
    test_gbk_genomes_mode_gene_test()
    test_gbk_genomes_mode_folder_fasta_only()
    test_gbk_genomes_mode_folder_fasta_only_cli()

