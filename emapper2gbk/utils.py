# Copyright (C) 2019-2024 Clémence Frioux & Arnaud Belcour - Inria Dyliss - Pleiade
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

import datetime
import itertools
import logging
import numpy as np
import os
import pandas as pd
import pronto
import requests
import simplejson
import shutil
import sys
import traceback

from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature as sf
from ete3 import NCBITaxa

try:
    # Import to be compatible with biopython version lesser than 1.78
    from Bio.Alphabet import IUPAC
except ImportError:
    # Exception to be compatible with biopython version superior to 1.78
    IUPAC = None

logger = logging.getLogger(__name__)


def get_basename(filepath):
    """Return the basename of given filepath.
    
    Args:
        filepath (str): path to a file
    
    Returns:
        str: basename

    >>> basename('~/an/interesting/file.txt')
    'file
    """
    return os.path.splitext(os.path.basename(filepath))[0]


def get_extension(filepath):
    """Get the extension of a filepath
    
    Args:
        filepath (str): path to a file
    
    Returns:
        str: extention of the file

    >>> extension('~/an/interesting/file.lp')
    'lp'
    >>> extension('nothing')
    ''
    >>> extension('nothing.important')
    'important'
    """
    return os.path.splitext(os.path.basename(filepath))[1][1:]


def is_valid_path(filepath):
    """Return True if filepath is valid.
    
    Args:
        filepath (str): path to file
    
    Returns:
        bool: True if path exists, False otherwise
    """
    if filepath and not os.access(filepath, os.W_OK):
        try:
            open(filepath, 'w').close()
            os.unlink(filepath)
            return True
        except OSError:
            return False
    else:  # path is accessible
        return True


def is_valid_file(filepath):
    """Return True if filepath exists.

    Args:
        filepath (str): path to file

    Returns:
        bool: True if path exists, False otherwise
    """
    try:
        open(filepath, 'r').close()
        return True
    except OSError:
        return False


def is_valid_dir(dirpath):
    """Return True if directory exists or can be created (then create it).
    
    Args:
        dirpath (str): path of directory

    Returns:
        bool: True if dir exists, False otherwise
    """
    if not os.path.isdir(dirpath):
        try:
            os.makedirs(dirpath)
            return True
        except OSError:
            return False
    else:
        return True

def check_valid_path(file_paths):
    """ Check validity of inputs.

    Args:
        file_paths (list): list of paths of input files to check

    """
    for elem in file_paths:
        if not is_valid_file(elem):
            logger.critical(f"{elem} is not a valid path file.")
            sys.exit(1)


def create_GO_namespaces_alternatives(gobasic_file=None):
    """
    Use pronto to query the Gene Ontology and to create the Ontology.
    Create a dictionary which contains for all GO terms their GO namespaces (molecular_function, ..).
    Create a second dictionary containing alternative ID for some GO terms (deprecated ones).

    Args:
        gobasic_file (str): path of GO basic file (if not provided will be downloaded)

    Returns:
        go_namespaces (dict): dictionary of GO terms namespace (key: GO Term ID, value: namespace associated to GO Term)
        go_alternatives (dict): dictionary of GO terms alternatives ID (key: GO Term ID, value: alternatives GO Term associated to GO Term)
    """
    go_basic_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'

    if gobasic_file:
        if not is_valid_file(gobasic_file):
            logger.critical(gobasic_file + ' file is not available, emapper2gbk will download the GO ontology and create it.')
            response = requests.get(go_basic_obo_url, stream=True)
            with open(gobasic_file, 'wb') as go_basic_file_write:
                shutil.copyfileobj(response.raw, go_basic_file_write)
    else:
        gobasic_file = go_basic_obo_url

    try:
        go_ontology = pronto.Ontology(gobasic_file)
    except SyntaxError as syntax_error:
        traceback.print_exc()
        logger.critical('/!\\ It seems that the obo fie (' + gobasic_file + ') has a syntax error.')
        logger.critical('You can try other obo file like the one from the github of gene Ontology (which can contains the fix for the syntax error): https://github.com/geneontology/go-ontology/blob/104252cc533a83248ed4df2caefd434e6eb8b531/src/ontology/go-edit.obo.')
        logger.critical('Or the obo file from the test folder of emapper2gbk (https://github.com/AuReMe/emapper2gbk/blob/master/tests/go-basic.obo).')
        sys.exit(1)

    # For each GO terms look to the namespaces associated with them.
    go_namespaces = {}
    for go_term in go_ontology:
        if 'GO:' in go_term:
            go = go_ontology[go_term]
            go_namespaces[go_term] = go.namespace

    # For each GO terms look if there is an alternative ID fo them.
    go_alternative = {}
    for go_term in go_ontology:
        if 'GO:' in go_term:
            go = go_ontology[go_term]
            if go.alternate_ids != frozenset():
                for go_alt in go.alternate_ids:
                    go_alternative[go_alt] = go_term

    return go_namespaces, go_alternative


def create_taxonomic_data(species_name):
    """
    Query the EBI with the species name to create a dictionary containing taxon id,
    taxonomy and some other informations.

    Args:
        species_name (str): species name (must be with genus for example "Escherichia coli")

    Returns:
        species_informations (dict): dictionary containing information about species
    """
    species_informations = {}

    compatible_species_name = species_name.replace('/', '_')

    if species_name == "bacteria":
        species_informations = {'db_xref': 'taxon:2', 'scientificName': 'Bacteria', 'commonName': 'eubacteria', 'formalName': 'false', 'rank': 'superkingdom', 'data_file_division': 'PRO', 'geneticCode': '11', 'submittable': 'false', 'description': 'bacteria genome', 'organism': 'bacteria', 'keywords': ['bacteria']} 
    elif species_name == "archaea":
        species_informations = {'db_xref': 'taxon:2157', 'scientificName': 'Archaea', 'formalName': 'false', 'rank': 'superkingdom', 'data_file_division': 'PRO', 'geneticCode': '11', 'submittable': 'false'}
    elif species_name == "eukaryota":
        species_informations = {'db_xref': 'taxon:2759', 'scientificName': 'Eukaryota', 'commonName': 'eucaryotes', 'formalName': 'false', 'rank': 'superkingdom', 'data_file_division': 'INV', 'geneticCode': '1', 'mitochondrialGeneticCode': '1', 'plastIdGeneticCode': '11', 'submittable': 'false'}
    elif species_name == "metagenome":
        species_informations = {'db_xref': 'taxon:256318', 'scientificName': 'metagenome', 'formalName': 'false', 'rank': 'species', 'division': 'UNC', 'lineage': 'unclassified sequences; metagenomes; ', 'geneticCode': '11', 'mitochondrialGeneticCode': '2', 'plastIdGeneticCode': '11', 'submittable': 'true'}
    elif species_name == 'cellular organisms':
        species_informations = {'db_xref': 'taxon:131567', 'scientificName': 'cellular organisms', 'formalName': 'false', 'rank': 'no rank', 'division': 'UNC', 'geneticCode': '1', 'submittable': 'false'}
    else:
        taxons = species_name.split(";")
        taxon_found = False

        for index, taxon in reversed(list(enumerate(taxons))):
            if not taxon or taxon == '':
                continue
            compatible_species_name = taxon.replace('/', '_')
            species_name_url = taxon.replace(' ', '%20')
            url = 'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/scientific-name/' + species_name_url

            try:
                response = requests.get(url)
            except requests.exceptions.ConnectionError:
                logger.critical('/!\\ No internet connection, check the connection or use the --ete option (if you have the NCBITaxa database already downloaded).')
                return None

            # Check if there is taxonomy information in the EBI response JSON.
            try:
                temp_species_informations = response.json()
                taxon_found = True
            except simplejson.errors.JSONDecodeError:
                if index==len(taxons):
                    logger.critical('/!\\ Error with {0} this taxa has not been found in https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/scientific-name/'.format(species_name))
                    logger.critical('/!\\ Check the name of the taxa and its presence in the EBI taxonomy database.')
                    logger.critical('/!\\ No genbank will be created for {0}.'.format(taxon))
                    return None
                else:
                    continue

            if temp_species_informations == []:
                logger.critical('/!\\ Error with {0} this taxa has not been found in https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/scientific-name/'.format(taxon))
                continue
            if taxon_found is True:
                break

        temp_species_informations = temp_species_informations[0]

        for temp_species_information in temp_species_informations:
            if temp_species_information == 'lineage':
                species_informations['taxonomy'] = temp_species_informations[temp_species_information].split('; ')[:-1]
            elif temp_species_information == 'division':
                species_informations['data_file_division'] = temp_species_informations[temp_species_information]
            elif temp_species_information == 'taxId':
                species_informations['db_xref'] = 'taxon:' + str(temp_species_informations[temp_species_information])
            else:
                species_informations[temp_species_information] = temp_species_informations[temp_species_information]

    species_informations['description'] = compatible_species_name + ' genome'
    species_informations['organism'] = compatible_species_name
    species_informations['keywords'] = [compatible_species_name]

    return species_informations


def create_taxonomic_data_ete(species_name):
    """
    Query ete taxonomy with the species name to create a dictionary containing taxon id,
    taxonomy and some other informations.
    Useful when no internet connection is available and the NCBITaxa database have already been downloaded.

    Args:
        species_name (str): species name (must be with genus for example "Escherichia coli")

    Returns:
        species_informations (dict): dictionary containing information about species
    """
    species_informations = {}
    taxons = species_name.split(";")
    for index, taxon in reversed(list(enumerate(taxons))):
        if not taxon:
            continue
        compatible_species_name = taxon.replace('/', '_')
        species_informations['description'] = compatible_species_name + ' genome'
        species_informations['organism'] = compatible_species_name
        species_informations['keywords'] = [compatible_species_name]

        ncbi = NCBITaxa()
        species_taxids = ncbi.get_name_translator([taxon])
        if taxon in species_taxids:
            species_taxid = species_taxids[taxon][-1]
            species_informations['db_xref'] = 'taxon:' + str(species_taxid)
            break
        else:
            if index ==0:
                logger.critical('/!\\ Error with {} this taxa has not been found in ete3 NCBITaxa Database'.format(taxon))
                logger.critical('/!\\ Check the name of the taxa and its presence in the NCBITaxa database.')
                logger.critical('/!\\ No genbank will be created for {}.'.format(taxon))
                return None
            else:
                continue

    return species_informations


def record_info(record_id, record_seq, species_informations):
    """ Create SeqBio record informations from species_informations dictionary and record id and record seq.

    Args:
        record_id (str): ID of the record (either contig/scaffold/chromosome ID)
        record_seq (Bio.Seq.Seq): Sequence associated to the record
        species_informations (dict): Dictionary containing information about species

    Returns:
        record (Bio.SeqRecord.SeqRecord): New record containing the informations provided
    """
    if record_id.isnumeric():
        newname = f"_{record_id}"
    elif "|" in record_id:
        newname = record_id.split("|")[1]
    else:
        newname = record_id

    record = SeqRecord(record_seq, id=record_id, name=newname,
                    description=species_informations['description'],
                    annotations={"molecule_type": "DNA"})

    # Condition to be compatible with biopython version lesser than 1.78
    if IUPAC:
        record.seq.alphabet = IUPAC.ambiguous_dna

    if 'data_file_division' in species_informations:
        record.annotations['data_file_division'] = species_informations['data_file_division']
    record.annotations['date'] = datetime.date.today().strftime('%d-%b-%Y').upper()
    if 'topology' in species_informations:
        record.annotations['topology'] = species_informations['topology']
    record.annotations['accessions'] = record_id
    if 'organism' in species_informations:
        record.annotations['organism'] = species_informations['organism']
    # Use of literal_eval for taxonomy and keywords to retrieve list.
    if 'taxonomy' in species_informations:
        record.annotations['taxonomy'] = species_informations['taxonomy']
    if 'keywords' in species_informations:
        record.annotations['keywords'] = species_informations['keywords']
    if 'source' in species_informations:
        record.annotations['source'] = species_informations['source']

    new_feature_source = sf.SeqFeature(sf.FeatureLocation(1-1,
                                                        len(record_seq)),
                                                        type="source")
    new_feature_source.qualifiers['scaffold'] = record_id
    if 'isolate' in species_informations:
        new_feature_source.qualifiers['isolate'] = species_informations['isolate']
    # db_xref corresponds to the taxon NCBI ID.
    # Important if you want to use Pathway Tools after.
    if 'db_xref' in species_informations:
        new_feature_source.qualifiers['db_xref'] = species_informations['db_xref']
    if 'cell_type' in species_informations:
        new_feature_source.qualifiers['cell_type'] = species_informations['cell_type']
    if 'dev_stage' in species_informations:
        new_feature_source.qualifiers['dev_stage'] = species_informations['dev_stage']
    if 'mol_type' in species_informations:
        new_feature_source.qualifiers['mol_type'] = species_informations['mol_type']

    record.features.append(new_feature_source)

    return record


def read_annotation(eggnog_outfile:str):
    """Read an eggnog-mapper annotation file and retrieve EC numbers and GO terms by genes.

    Args:
        eggnog_outfile (str): path to eggnog-mapper annotation file

    Returns:
        dict: dict of genes and their annotations as {gene1:{EC:'..,..', GOs:'..,..,'}}
    """
    # Look at the twentieth first rows to find the header.
    with open(eggnog_outfile, 'r') as f:
        twentieth_first_rows = list(itertools.islice(f, 20))
        first_row_after_header = min([index for index, str_row in enumerate(twentieth_first_rows) if not str_row.startswith('#')])
        header_row = first_row_after_header - 1
        headers_row = twentieth_first_rows[header_row].lstrip("#").strip().split('\t')

    # Fix issue when header is incomplete (eggnog before version 2.0).
    if len(headers_row) == 17:
        headers_row.extend(['tax_scope', 'eggNOG_OGs', 'bestOG', 'COG_functional_category', 'eggNOG_free_text'])

    to_extract_annotations = ['GOs','EC', 'Preferred_name']
    if 'PFAMs' in headers_row:
        to_extract_annotations.append('PFAMs')
    if 'BiGG_Reaction' in headers_row:
        to_extract_annotations.append('BiGG_Reaction')
    if 'KEGG_Reaction' in headers_row:
        to_extract_annotations.append('KEGG_Reaction')
    if 'CAZy' in headers_row:
        to_extract_annotations.append('CAZy')

    # Use chunk when reading eggnog file to cope with big file.
    chunksize = 10 ** 6
    for annotation_data in pd.read_csv(eggnog_outfile, sep='\t', comment='#', header=None, dtype = str, chunksize = chunksize):
        annotation_data.replace(np.nan, '', inplace=True)
        # Assign the headers
        annotation_data.columns = headers_row
        if 'query_name' in annotation_data.columns:
            # Check if the gene IDs are numeric, if yes add 'gene_' in front of them.
            numeric_row_dataframe = pd.to_numeric(annotation_data['query_name'], errors='coerce').notnull()
            if bool(numeric_row_dataframe.any()) is True:
                annotation_data.loc[numeric_row_dataframe, 'query_name'] = 'gene_' + annotation_data.loc[numeric_row_dataframe, 'query_name']
            annotation_dict = annotation_data.set_index('query_name')[to_extract_annotations].to_dict('index')
        # 'query' added for compatibility with eggnog-mapper 2.1.2
        elif 'query' in annotation_data.columns:
            numeric_row_dataframe = pd.to_numeric(annotation_data['query'], errors='coerce').notnull()
            if bool(numeric_row_dataframe.any()) is True:
                annotation_data.loc[numeric_row_dataframe, 'query'] = 'gene_' + annotation_data.loc[numeric_row_dataframe, 'query']
            annotation_dict = annotation_data.set_index('query')[to_extract_annotations].to_dict('index')
        for key in annotation_dict:
            yield key, annotation_dict[key]


def create_cds_feature(id_gene, start_position, end_position,
                        strand, annot, go_namespaces, go_alternatives,
                        gene_protein_seq, gff_extracted_annotations=None,
                        location_exons=None):
    """ Create Biopython CDS feature from gene ID, gene positions, gene sequecne and gene annotations.

    Args:
        id_gene (str): ID of the gene
        start_position (int): start position of the gene
        end_position (int): end position of the gene
        strand (str): strand of teh gene
        annot (dict): dictionary of eggnog-ammper annotation (key: gene_id, value: ['GOs','EC', 'Preferred_name'])
        go_namespaces (dict): dictionary of GO terms namespace (key: GO Term ID, value: namespace associated to GO Term)
        go_alternatives (dict): dictionary of GO terms alternatives ID (key: GO Term ID, value: alternatives GO Term associated to GO Term)
        gene_protein_seq (dict): dictionary of protein sequence associated to genes (key: gene id, value: sequence)
        gff_extracted_annotations (dict): dictionary of annotation to copy from the GFF file into the Genbank file.
        location_exons (list): list of FeatureLocation containing exons positions for the CDS, None if no exon position or if not in gff type gmove.

    Returns:
        new_feature_cds (Bio.SeqFeature.SeqFeature): New SeqFeature containing the informations provided
    """
    if location_exons and len(location_exons)>=2:
        exon_feature_locations = sf.CompoundLocation(location_exons, operator='join')

        new_feature_cds = sf.SeqFeature(exon_feature_locations, type='CDS')
    else:
        new_feature_cds = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                            end_position,
                                                            strand),
                                                                    type="CDS")

    new_feature_cds.qualifiers['locus_tag'] = id_gene

    if gff_extracted_annotations:
        if 'product' in gff_extracted_annotations:
            gff_product = gff_extracted_annotations['product']
        else:
            gff_product = None
    else:
        gff_product = None

    # Add GO annotation according to the namespace.
    if id_gene in annot.keys():
        # Add gene name.
        if 'Preferred_name' in annot[id_gene]:
            if annot[id_gene]['Preferred_name'] != '' and annot[id_gene]['Preferred_name'] != '-':
                new_feature_cds.qualifiers['gene'] = annot[id_gene]['Preferred_name']

        # Add product name from GFF.
        if gff_product:
            new_feature_cds.qualifiers['product'] = gff_product

        # Add GO terms.
        if 'GOs' in annot[id_gene]:
            gene_gos = annot[id_gene]['GOs'].split(',')
            if '' in gene_gos:
                gene_gos.remove('')
            # '-' added for compatibility with eggnog-mapper 2.1.2
            if '-' in gene_gos:
                gene_gos.remove('-')
            if gene_gos != []:
                go_components = []
                go_functions = []
                go_process = []

                for go in gene_gos:
                    # Check if GO term is not a deprecated one.
                    # If yes take the corresponding one in alternative GO.
                    if go not in go_namespaces:
                        go_test = go_alternatives[go]
                    else:
                        go_test = go
                    if go_namespaces[go_test] == 'cellular_component':
                            go_components.append(go)
                    if go_namespaces[go_test] == 'molecular_function':
                        go_functions.append(go)
                    if go_namespaces[go_test] == 'biological_process':
                        go_process.append(go)
                new_feature_cds.qualifiers['go_component'] = go_components
                new_feature_cds.qualifiers['go_function'] = go_functions
                new_feature_cds.qualifiers['go_process'] = go_process


        # Add EC annotation.
        if 'EC' in annot[id_gene]:
            gene_ecs = annot[id_gene]['EC'].split(',')
            if '' in gene_ecs:
                gene_ecs.remove('')
            if '-' in gene_ecs:
                gene_ecs.remove('-')
            if gene_ecs != []:
                new_feature_cds.qualifiers['EC_number'] = gene_ecs

        # Add Pfam dbxref.
        if 'PFAMs' in annot[id_gene]:
            gene_pfams = annot[id_gene]['PFAMs'].split(',')
            if '' in gene_pfams:
                gene_pfams.remove('')
            if '-' in gene_pfams:
                gene_pfams.remove('-')
            if gene_pfams != []:
                new_feature_cds.qualifiers['dbxref'] = ['PFAM:'+pfam for pfam in gene_pfams]

        # Add CAZ dbxref.
        if 'CAZy' in annot[id_gene]:
            gene_cazys = annot[id_gene]['CAZy'].split(',')
            if '' in gene_cazys:
                gene_cazys.remove('')
            if '-' in gene_cazys:
                gene_cazys.remove('-')
            if gene_cazys != []:
                new_feature_cds.qualifiers['dbxref'] = ['CAZY:'+kegg for kegg in gene_cazys]

        # Add bigg dbxref.
        if 'BiGG_Reaction' in annot[id_gene]:
            gene_reaction_biggs = annot[id_gene]['BiGG_Reaction'].split(',')
            if '' in gene_reaction_biggs:
                gene_reaction_biggs.remove('')
            if '-' in gene_reaction_biggs:
                gene_reaction_biggs.remove('-')
            if gene_reaction_biggs != []:
                new_feature_cds.qualifiers['dbxref'] = ['BIGG:'+bigg for bigg in gene_reaction_biggs]

        # Add kegg dbxref.
        if 'KEGG_Reaction' in annot[id_gene]:
            gene_reaction_keggs = annot[id_gene]['KEGG_Reaction'].split(',')
            if '' in gene_reaction_keggs:
                gene_reaction_keggs.remove('')
            if '-' in gene_reaction_keggs:
                gene_reaction_keggs.remove('-')
            if gene_reaction_keggs != []:
                new_feature_cds.qualifiers['dbxref'] = ['KEGG:'+kegg for kegg in gene_reaction_keggs]

    if id_gene in gene_protein_seq:
        new_feature_cds.qualifiers['translation'] = gene_protein_seq[id_gene]

    return new_feature_cds
