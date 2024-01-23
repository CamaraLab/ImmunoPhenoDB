import requests
import json
import logging
import pandas as pd
import csv
import urllib
from datetime import datetime
import warnings

import configparser
import mysql.connector # mysql-connector-python
from mysql.connector import errorcode

import sys
from tqdm.autonotebook import tqdm
from unittest.mock import patch
from contextlib import contextmanager

SCI_CRUNCH_BASE = "http://www.scicrunch.org"
SCI_RRID_ENDPOINT = "/resolver/"
SCI_FILE = ".json"
UNIPROT_BASE = "https://rest.uniprot.org"
UNIPROT_ENDPOINT = "/uniprotkb/search"
EBI_BASE = "https://www.ebi.ac.uk/ols4/api/search"

def _read_cells(csv_file: str) -> list:
    """
    Converts a csv to a pandas dataframe and returns a list of cell names

    Parameters:
        csv_file (str): name of csv file containing normalized counts

    Returns:
        list containing all cell names
    """
    normalized_counts = pd.read_csv(csv_file, header=0, sep=",", index_col=[0])
    return list(normalized_counts.index)

def _read_experiment(csv_file: str) -> dict:
    """ 
    Parses a spreadsheet containing information about an experiment

    Parameters:
        csv_file (str): name of csv file containing experiment details
    
    Returns:
        experiment (dict): dictionary containing information about the
            experiment's name, type, pmid, doi, and tissue type
    """
    experiment = {}
    expTypes = ['cite', 'reap', 'ab']
 
    with open(csv_file, 'r') as csv_file:
        reader = csv.reader(csv_file)

        for row in reader:
            if 'name experiment' in row[0].lower():
                experiment['exp'] = row[1]
            
            if 'type' in row[0].lower():
                if any(expType in row[1].lower() for expType in expTypes):
                    experiment['type'] = row[1]
            
            if 'pmid' in row[0].lower():
                experiment['pmid'] = row[1]

            if 'doi' in row[0].lower():
                experiment['doi'] = row[1]
            
            if 'tissue' in row[0].lower():
                experiment['tissue'] = row[1]

            if 'antibody table' in row[0].lower():
                break
                
    return experiment

def _read_antibodies(csv_file: str) -> list:
    """
    Parses a spreadsheet containing antibody IDs 

    Parameters:
        csv_file (str): name of csv file containing antibody names
    
    Returns:
        antibodies (list): list of antibody ids (ex: AB_XXXXXXX)
    """
    antibodies = []

    with open(csv_file, 'r') as csv_file:
        lines = csv_file.readlines()
        num_lines = len(lines)
        ab_index = (lines.index('Antibody table:,\n') + 1)
        
        while ab_index < num_lines:
            ab = lines[ab_index].strip().split(",", maxsplit=1)
            ab_name_strip = ab[0].strip()
            ab_id_strip = ab[1].strip()
            antibodies.append([ab_name_strip, ab_id_strip])
            ab_index += 1
    
    return antibodies

def filter_antibodies(protein_filepath: pd.DataFrame,
                      csv_file: str) -> pd.DataFrame:
    """
    Filters ADT protein table using only antibodies found in 
    a user-provided spreadsheet
    
    Parameters:
        protein_filepath (str): file path to protein data
        csv_file (str): file path to provided spreadsheet
    
    Returns:
        filt_df (pd.DataFrame): dataframe containing rows
            that reflect antibodies listed in the spreadsheet
    """
    
    antibody_pairs = _read_antibodies(csv_file)
    antibodies_list = [ab[0] for ab in antibody_pairs]

    protein_data = pd.read_csv(protein_filepath, sep=',', index_col=[0])
    
    # Subset the columns that are in our spreadsheet
    filt_df = protein_data.loc[antibodies_list]
    
    return filt_df

def _ab_dict(specs_csv: str) -> dict:
    """
    Constructs a dictionary of each antibody name with their id

    Parameters:
        specs_csv (str): name of csv file containing antibodies and ids

    Returns:
        ab_dict (dict): dictionary containing antibody names as keys, ids as values
    """
    antibodies = _read_antibodies(specs_csv)

    ab_dict = {}

    for ab in antibodies:
        ab_dict[ab[0]] = ab[1]

    return ab_dict

def _findkeys(node: list, kv: str) -> iter:
    """
    Recursively find all instances of a key in a nested JSON file

    Parameters:
        node (list/dict): root node (key) in JSON file or dictionary
        kv (str): search term for key

    Return:
        x (iterator): all key-value results
    """
    if isinstance(node, list):
        for i in node:
            for x in _findkeys(i, kv):
                yield x
    elif isinstance(node, dict):
        if kv in node:
            yield node[kv]
        for j in node.values():
            for x in _findkeys(j, kv):
                yield x

@contextmanager
def _tqdm_output(tqdm, write=sys.stderr.write):
    """
    Context manager wrapper for tqdm progress bar

    Parameters:
        tqdm (tqdm object): object created with the tqdm class
        write (file output object): location to display progress bar
    
    Returns:
        tqdm iterator when using a context manager
    """
    def wrapper(message):
        if message != '\n':
            tqdm.clear()
        write(message)
        if '\n' in message:
            tqdm.display()

    with patch('sys.stdout', sys.stderr), patch('sys.stderr.write', wrapper):
        yield tqdm

def _uniprot(sci_crunch_alias: str = None, 
             user_uniprotID: str = None) -> tuple:
    """
    Retrieves information about an antibody using UniProt's API

    Parameters:
        sci_crunch_alias (str): antibody alias returned from SciCrunch API
        uniprotID (str): user-provided UniProtID for manual entry

    Returns:
        uniprotID, otherAliases (tuple [str, list]): antibody alias
            and aliases retrieved from UniProt
    """

    otherAliases = []

    # Protein and gene search (strictest level)
    protein_gene_params = {
        'query':  f'protein_name:{sci_crunch_alias} AND gene_exact:{sci_crunch_alias} AND organism_id:9606 AND reviewed:true',
        'fields': 'gene_primary, gene_synonym, protein_name',
        'format': 'json'
    }

    # Protein search query params
    protein_params = {
        'query':  f'protein_name:{sci_crunch_alias} AND organism_id:9606 AND reviewed:true',
        'fields': 'gene_primary, gene_synonym, protein_name',
        'format': 'json'
    }

    # Gene search query params if protein query fails (response != 200 or no results)
    gene_params = {
        'query':  f'gene_exact:{sci_crunch_alias} AND organism_id:9606 AND reviewed:true',
        'fields': 'gene_primary, gene_synonym, protein_name',
        'format': 'json'
    }

    # Direct UniProtID accession search query params if user provided UniProtID
    accession_params = {
        'query':  f'accession:{user_uniprotID}',
        'fields': 'gene_primary, gene_synonym, protein_name',
        'format': 'json'
    }

    # Searching for UniProtID using SciCrunch alias if not user provided
    if sci_crunch_alias is not None and user_uniprotID is None:
        # Strictest level search
        protein_gene_response = requests.get(UNIPROT_BASE + UNIPROT_ENDPOINT, params=protein_gene_params)
        protein_gene_JSON = protein_gene_response.json()
        
        if protein_gene_response.status_code == 200 and len(protein_gene_JSON['results']) != 0:
            resJSON = protein_gene_JSON
        else:
            # Protein only search
            protein_response = requests.get(UNIPROT_BASE + UNIPROT_ENDPOINT, params=protein_params)
            protein_JSON = protein_response.json()

            if protein_response.status_code == 200 and len(protein_JSON['results']) != 0:
                resJSON = protein_JSON
            else:
                # Gene only search
                gene_response = requests.get(UNIPROT_BASE + UNIPROT_ENDPOINT, params=gene_params)
                gene_JSON = gene_response.json()
                resJSON = gene_JSON

    # Searching for UniProtID if accession ID is user provided
    elif sci_crunch_alias is None and user_uniprotID is not None:
        accession_response = requests.get(UNIPROT_BASE + UNIPROT_ENDPOINT, params=accession_params)
        resJSON = accession_response.json()
        
    # Retrieve Uniprot ID using SciCrunch alias (required)
    try:
        uniprotID = resJSON['results'][0]['primaryAccession']
    except:
        raise Exception(f"Unable to find UniProtID for: {sci_crunch_alias}")

    # Retrieve recommended protein name (required)
    try:
        recommendedName = (resJSON['results'][0]
                            ['proteinDescription']
                            ['recommendedName']
                            ['fullName']
                            ['value'])
        otherAliases.append(recommendedName)
    except:
        raise Exception(f"Unable to find UniProt recommended name for: {sci_crunch_alias}")
    
    # Retrieve alternative protein names (optional)
    try:
        alternativeNamesDicts = (resJSON['results'][0]
                            ['proteinDescription']
                            ['alternativeNames'])
        alternativeNames = list(_findkeys(alternativeNamesDicts, 'value'))
        otherAliases.extend(alternativeNames)
    except:
        pass

    # Retrieve CD antigen names (optional)
    try:
        # CD Antigen names
        antigenNamesDicts = (resJSON['results'][0]
                            ['proteinDescription']
                            ['cdAntigenNames'])
        antigenNames = list(_findkeys(antigenNamesDicts, 'value'))
        otherAliases.extend(antigenNames)
    except:
        pass

    # Retrieve gene names (optional)
    try:
        # Gene names
        geneName = (resJSON['results'][0]
                            ['genes'][0]
                            ['geneName']
                            ['value'])
        otherAliases.append(geneName)
    except:
        pass

    # Retrieve gene name synonyms (optional)
    try:
        geneSynonymsDict = (resJSON['results'][0]
                            ['genes'][0]
                            ['synonyms'])
        geneSynynoms = list(_findkeys(geneSynonymsDict, 'value'))
        otherAliases.extend(geneSynynoms)      
    except:
        pass

    return uniprotID, otherAliases

def _sci_crunch_hits(ab_id: str) -> bool:
    """
    Retrieves information about an antibody using SciCrunch's API

    Parameters:
        ab_id (str): antibody id from The Antibody Registry
    
    Returns:
        response (json): JSON result if successful response
        false (bool): False if error with response
    """
    response = requests.get(SCI_CRUNCH_BASE + SCI_RRID_ENDPOINT + ab_id + SCI_FILE)
    if response.status_code == 200:
        return response.json()
    else:
        return False

def _sci_uni(ab_id_pair: list, 
             user_uniprotID: str = None) -> list:
    """
    Retrieves information from SciCrunch and UniProt for an antibody

    Parameters:
        ab_id_pair (list[str, str]): antibody name and id ('CD90', 'AB_123')
        user_uniprotID (str): UniProt accession ID (optional, for manual entry)

    Returns:
        each_hit_results (list): information gathered using SciCrunch and UniProt
            for an antibody.
            Example: [(alias, polyclonal, host, uniprotID, otherAliases), ...]
    """
    errors = []

    each_hit_results = []

    csv_alias = ab_id_pair[0]
    csv_ab_id = ab_id_pair[1]

    ### SciCrunch ###
    sciJSON = _sci_crunch_hits(csv_ab_id)

    if sciJSON != False:
        num_hits = sciJSON['hits']['total']

        for i in range(num_hits):
            sci_alias = (sciJSON['hits']['hits'][i]
                            ['_source']['antibodies']['primary'][0]
                            ['targets'][0]
                            ['name'])

            clonality = (sciJSON['hits']['hits'][i]
                                ['_source']['antibodies']['primary'][0]
                                ['clonality']['name'])
            
            if clonality != 'polyclonal':
                clonalBool = False

            cloneID = (sciJSON['hits']['hits'][i]
                            ['_source']['antibodies']['primary'][0]
                            ['clone']['identifier']).replace('Clone ', '')
            
            # If there's a cloneID for polyclonal, raise an error
            # Reverse: If monoclonal, and it does NOT have a clone ID: then spit an error
            if clonalBool and cloneID:
                errors.append((sci_alias, csv_ab_id, 'CloneID found for polyclonal antibody'))
                logging.warning(f"Skipping {sci_alias}. Refer to error log.")
                continue # skip this antibody hit

            host = (sciJSON['hits']['hits'][i]
                        ['_source']['organisms']['source'][0]
                        ['species']['name'])

            ### UniProt ###
            try:
                if user_uniprotID is None:
                    uniprotID, otherAliases = _uniprot(sci_crunch_alias=sci_alias)
                    each_hit_results.append((sci_alias, clonalBool, host, uniprotID, cloneID, otherAliases))
                elif user_uniprotID is not None:
                    uniprotID, otherAliases = _uniprot(user_uniprotID=user_uniprotID)
                    each_hit_results.append((sci_alias, clonalBool, host, uniprotID, cloneID, otherAliases))
            except Exception as e:
                errors.append((ab_id_pair[0], ab_id_pair[1], str(e)))
                logging.warning(f"Skipping {csv_alias}. Refer to antibody_errors.txt")
                continue # skip this antibody hit
    else:
        errors.append((ab_id_pair[0], ab_id_pair[1], 
                        f'Unable to find antibody ID: {ab_id_pair[1]} in SciCrunch'))
        logging.warning(f"Skipping {ab_id_pair[0]}. Refer to antibody_errors.txt")
 
    # Log errors to an external file
    if errors:
        with open('antibody_errors.txt', 'a') as f:
            for error in errors:
                time = str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
                f.write(f"[{time}] {error}\n")

    return each_hit_results

def _config(filename: str = 'config.ini', 
            section: str ='mysql') -> dict:
    """
    Configures authorization parameters for a MySQL Database

    Parameters:
        filename (str): .ini file containing: host, user, password, database
        section (str): type of database "[mysql]" found at beginning of .ini file
    
    Returns:
        db (dict): dictionary containing login credentials such as host, 
        user, password, database, and port
    """
    parser = configparser.ConfigParser()
    parser.read(filename)

    db = {}

    if parser.has_section(section):
        if parser.has_section(section):
            params = parser.items(section)

        for param in params:
            db[param[0]] = param[1]
    else:
        raise Exception(f"{section} not found in {filename}")
  
    return db

def _connect_db_tables():
    """
    SQL Queries used to generate a schema containing the following tables:
        proteins
        aliases
        antibodies
        antigens
        tissues
        cell_types
        experiments
        cells
        antigen_expression
    """
    proteins_table = """CREATE TABLE IF NOT EXISTS proteins(
                    idUniProtKB CHAR(6) NOT NULL,
                    PRIMARY KEY (idUniProtKB));
                    """

    aliases_table = """CREATE TABLE IF NOT EXISTS aliases(
                    alias CHAR(255) NOT NULL,
                    idUniProtKB CHAR(6) NOT NULL,
                    PRIMARY KEY (alias, idUniProtKB),
                    FOREIGN KEY (idUniProtKB) REFERENCES proteins(idUniProtKB));
                    """

    antibodies_table = """CREATE TABLE IF NOT EXISTS antibodies(
                        idAntibody CHAR(11) NOT NULL,
                        polyclonal BOOLEAN NOT NULL,
                        host VARCHAR(45) NOT NULL,
                        cloneID VARCHAR(45) NOT NULL,
                        abName VARCHAR(128) NOT NULL,
                        abTarget VARCHAR(128) NOT NULL,
                        citation VARCHAR(128),
                        comments VARCHAR(2000),
                        vendor VARCHAR(128),
                        catalogNum VARCHAR(128),
                        FULLTEXT (idAntibody, host, cloneID, abName, abTarget, citation, comments, vendor, catalogNum),
                        PRIMARY KEY (idAntibody));
                        """

    antigens_table = """CREATE TABLE IF NOT EXISTS antigens(
                    idAntibody CHAR(11) NOT NULL,
                    idUniProtKB CHAR(6) NOT NULL,
                    PRIMARY KEY (idAntibody, idUniProtKB),
                    FOREIGN KEY (idAntibody) REFERENCES antibodies(idAntibody),
                    FOREIGN KEY (idUniProtKB) REFERENCES proteins(idUniProtKB));
                    """
    
    tissues_table = """CREATE TABLE IF NOT EXISTS tissues (
                    idBTO CHAR(11) NOT NULL,
                    tissueName VARCHAR(128) NOT NULL,
                    PRIMARY KEY (idBTO));
                    """
    
    cell_types_table = """CREATE TABLE IF NOT EXISTS cell_types (
                            idCL CHAR(10) NOT NULL,
                            label VARCHAR(128) NOT NULL,
                            FULLTEXT (idCL, label),
                            PRIMARY KEY (idCL));
                            """
    
    experiments_table = """CREATE TABLE IF NOT EXISTS experiments(
                        idExperiment INT UNSIGNED NOT NULL AUTO_INCREMENT,
                        nameExp CHAR(255) NOT NULL,
                        typeExp ENUM('CITE', 'REAP', 'AB') NOT NULL,
                        pmid INT UNSIGNED NULL,
                        doi CHAR(255) NULL,
                        idBTO CHAR(11) NOT NULL,
                        PRIMARY KEY (idExperiment),
                        FOREIGN KEY (idBTO) REFERENCES tissues(idBTO));
                        """
    
    cells_table = """CREATE TABLE IF NOT EXISTS cells(
                    idCL CHAR(10) NULL, 
                    idCell INT UNSIGNED NOT NULL AUTO_INCREMENT,
                    idCellOriginal VARCHAR(255) NOT NULL,
                    idExperiment INT UNSIGNED NOT NULL,
                    certainty FLOAT NULL,
                    PRIMARY KEY (idCell, idExperiment),
                    FOREIGN KEY (idExperiment) REFERENCES experiments(idExperiment),
                    FOREIGN KEY (idCL) REFERENCES cell_types(idCL));
                    """
    
    antigen_expr_table = """CREATE TABLE IF NOT EXISTS antigen_expression(
                        idCell INT UNSIGNED NOT NULL, 
                        idAntibody CHAR(11) NOT NULL,
                        idExperiment INT UNSIGNED NOT NULL,
                        rawValue FLOAT NOT NULL,
                        normValue FLOAT NOT NULL,
                        background BOOL NOT NULL,
                        PRIMARY KEY (idCell, idAntibody, idExperiment),
                        FOREIGN KEY (idCell) REFERENCES cells(idCell),
                        FOREIGN KEY (idAntibody) REFERENCES antibodies(idAntibody),
                        FOREIGN KEY (idExperiment) REFERENCES cells(idExperiment));
                        """

    check_tables = """SHOW TABLES;"""

    # Connect to database
    params = _config()
    print('Connecting to the MySQL database...')
    conn = mysql.connector.connect(**params)
    print("Connected to db\n")
    cursor = conn.cursor()

    # Check for existing tables
    cursor.execute(check_tables)
    current_tables = cursor.fetchall()
    tables_list = [table[0] for table in current_tables]
    
    required_tables = ['aliases', 'antibodies', 'antigen_expression', 'antigens',
                       'cells', 'experiments', 'proteins', 'tissues', 'cell_types']

    if set(required_tables) == set(tables_list):
        print("All tables present...")
    else:
        print("Creating proteins table...")
        print("Creating aliases table...")
        print("Creating antibodies table...")
        print("Creating antigens table...")
        print("Creating tissues table...")
        print("Creating cell_types table...")
        print("Creating experiments table...")
        print("Creating cells table...")
        print("Creating antigen expression table...")

    # Create proteins table
    cursor.execute(proteins_table)
    conn.commit()

    # Create aliases table
    cursor.execute(aliases_table)
    conn.commit()

    # Create antibodies table
    cursor.execute(antibodies_table)
    conn.commit()

    # Create antigens table
    cursor.execute(antigens_table)
    conn.commit()

    # Create tissues table
    cursor.execute(tissues_table)
    conn.commit()

    # Create cell_types table
    cursor.execute(cell_types_table)
    conn.commit()

    # Create experiments table
    cursor.execute(experiments_table)
    conn.commit()

    # Create cells table
    cursor.execute(cells_table)
    conn.commit()

    # Create antigen-expression table
    cursor.execute(antigen_expr_table)
    conn.commit()

    if conn is not None:
        print("Disconnecting from database...\n")
        cursor.close()
        conn.close()

def _connect_db_procedures():
    """
    SQL Queries to generate MySQL procedures for inserting data into tables.
    Contains the following procedures:
        insert_experiment
        insert_antibody
        insert_alias
        insert_cell
        insert_antigen_expression
    """
    delete_experiment_proc = "DROP PROCEDURE IF EXISTS insert_experiment;"
    insert_experiment_proc = """CREATE PROCEDURE insert_experiment(
            IN nameExp CHAR(255),
            IN typeExp ENUM('CITE', 'REAP', 'AB'),
            IN pmid INT,
            IN doi CHAR(255),
            IN idBTO CHAR(11),
            IN tissueName VARCHAR(128) 
        ) 
        BEGIN
            INSERT IGNORE INTO tissues(idBTO, tissueName) VALUES (idBTO, tissueName);
            INSERT IGNORE INTO experiments(nameExp, typeExp, pmid, doi, idBTO) VALUES (nameExp, typeExp, pmid, doi, idBTO);
        END
    """

    delete_antibody_proc = """DROP PROCEDURE IF EXISTS insert_antibody;"""
    insert_antibody_proc = """CREATE PROCEDURE insert_antibody(
            IN idAntibody CHAR(11),
            IN alias CHAR(255),
            IN idUniProtKB CHAR(6),
            IN polyclonal BOOLEAN,
            IN host VARCHAR(45),
            IN cloneID VARCHAR(45),
            IN abName VARCHAR(128),
            IN abTarget VARCHAR(128),
            IN citation VARCHAR(128),
            IN comments VARCHAR(2000),
            IN vendor VARCHAR(128),
            IN catalogNum VARCHAR(128)
        )
        BEGIN
            INSERT IGNORE INTO antibodies (idAntibody, polyclonal, host, cloneID, abName, abTarget, citation, comments, vendor, catalogNum) 
            VALUES (idAntibody, polyclonal, host, cloneID, abName, abTarget, citation, comments, vendor, catalogNum);
            
            INSERT IGNORE INTO proteins (idUniProtKB) VALUES (idUniProtKB);
            INSERT IGNORE INTO aliases (alias, idUniProtKB) VALUES (alias, idUniProtKB);
            INSERT IGNORE INTO antigens (idAntibody, idUniProtKB) VALUES (idAntibody, idUniProtKB);
        END
    """

    delete_alias_proc = """DROP PROCEDURE IF EXISTS insert_alias;"""
    insert_alias_proc = """CREATE PROCEDURE insert_alias(
            IN alias CHAR(255),
            IN idUniProtKB CHAR(6)
        )
        BEGIN
            INSERT IGNORE INTO aliases(alias, idUniProtKB) VALUES (alias, idUniProtKB);
        END 
    """

    delete_cell_proc = """DROP PROCEDURE IF EXISTS insert_cell;"""
    insert_cell_proc = """CREATE PROCEDURE insert_cell(
            IN idCL CHAR(10),
            IN idCellOriginal VARCHAR(255),
            IN idExperiment INT,
            IN certainty FLOAT,
            IN label VARCHAR(128)
        )
        BEGIN 
            INSERT IGNORE INTO cell_types(idCL, label) VALUES (idCL, label);
        
            INSERT INTO cells(idCL, idCellOriginal, idExperiment, certainty) 
            VALUES (idCL, idCellOriginal, idExperiment, certainty)
            ON DUPLICATE KEY UPDATE idCL=idCL, idCellOriginal=idCellOriginal, idExperiment=idExperiment, certainty=certainty;
        END 
    """

    delete_antigen_proc = """DROP PROCEDURE IF EXISTS insert_antigen_expression;"""
    insert_antigen_proc = """CREATE PROCEDURE insert_antigen_expression(
            IN idCell INT,
            IN idAntibody CHAR(11),
            IN idExperiment INT UNSIGNED,
            IN rawValue FLOAT,
            IN normValue FLOAT,
            IN background BOOLEAN
        )
        BEGIN
            INSERT IGNORE INTO antigen_expression(idCell, idAntibody, idExperiment, rawValue, normValue, background)
            VALUES (idCell, idAntibody, idExperiment, rawValue, normValue, background);
        END
    """
    check_proc = """SELECT specific_name
                    FROM information_schema.ROUTINES
                    WHERE routine_schema=DATABASE();"""
    
    # Connect to database
    params = _config()
    print('Connecting to the MySQL database...')
    conn = mysql.connector.connect(**params)
    print("Connected to db\n")
    cursor = conn.cursor()

    # Check for existing tables
    cursor.execute(check_proc)
    current_procs = cursor.fetchall()
    procs_list = [proc[0] for proc in current_procs]

    required_procs = ['insert_alias', 'insert_antibody', 'insert_antigen_expression',
                      'insert_cell', 'insert_experiment']

    if set(required_procs) == set(procs_list):
        print("All procedures present...")
    else:
        print("Creating insert_experiment procedure...")
        print("Creating insert_antibody procedure...")
        print("Creating insert_alias procedure...")
        print("Creating insert_cell procedure...")
        print("Creating insert_antigen procedure...")

    cursor.execute(delete_experiment_proc)
    conn.commit()

    cursor.execute(insert_experiment_proc)
    conn.commit()

    cursor.execute(delete_antibody_proc)
    conn.commit()

    cursor.execute(insert_antibody_proc)
    conn.commit()

    cursor.execute(delete_alias_proc)
    conn.commit()

    cursor.execute(insert_alias_proc)
    conn.commit()

    cursor.execute(delete_cell_proc)
    conn.commit()

    cursor.execute(insert_cell_proc)
    conn.commit()

    cursor.execute(delete_antigen_proc)
    conn.commit()
    
    cursor.execute(insert_antigen_proc)
    conn.commit()

    if conn is not None:
        print("Disconnecting from database...\n")
        cursor.close()
        conn.close()

def extra_ab_info(ab_id: str) -> pd.DataFrame:
    """
    Retrieves additional information about an antibody using SciCrunch

    Parameters:
        ab_id (str): antibody id
    
    Returns:
        final_df (pd.DataFrame): dataframe containing various
            information about an antibody (target, clonality, etc)
    """
    encoding = urllib.parse.quote(ab_id, safe='', encoding=None, errors=None)
    url = f"https://scicrunch.org/php/data-federation-csv.php?orMultiFacets=true&count=10000&nifid=nif-0000-07730-1&exportType=data&q={encoding}&offset=0"
    
    raw_df = pd.read_csv(url)
    
    # Select columns to add to resulting dataframe
    cols = ['id', 
            'ab_name', 
            'ab_target', 
            'proper_citation', 
            'clonality', 
            'defining_citation', 
            'comments', 
            'clone_id', 
            'target_species', 
            'vendor', 
            'catalog_num']
    
    final_df = raw_df[cols]
    
    return final_df

def _connect_db_antibody(specs_csv: str, IPD):
    """
    Connects to a MySQL database and inserts information about antibodies

    Parameters:
        specs_csv (str): name of csv file containing a spreadsheet with
            information about the experiment and antibodies
    """
    # Connect to database
    params = _config()
    print('Connecting to the MySQL database...')
    conn = mysql.connector.connect(**params)
    print("Connected to db\n")
    cursor = conn.cursor()

    print("Inserting antibodies...")

    errors = []
    
    # Read all antibodies from csv file ('CD90', 'AB_123')
    antibodies = _read_antibodies(specs_csv)

    with _tqdm_output(tqdm(antibodies, total=len(IPD.protein.columns))) as tqdm_ab:
        for ab_id_pair in tqdm_ab:
            # CHECK: Does this antibody id exist in the database?
            check_ab_exists_query = """SELECT COUNT(*)
                                    FROM antibodies 
                                    WHERE antibodies.idAntibody=(%s)"""
            # Parameter must be converted from str to tuple
            cursor.execute(check_ab_exists_query, (ab_id_pair[1], ))
            ab_exists_result = cursor.fetchone()[0]

            if ab_exists_result == 1 or ab_id_pair[0] not in IPD._protein_matrix.columns:
                # Antibody already exists, we can skip 
                continue
            else:
                # First check if this antibody is completely background       
                all_bg = all(x == 0 for x in (IPD.classified_filt.loc[:, ab_id_pair[0]]))
                if all_bg:
                    # Skip this antibody
                    errors.append((ab_id_pair[0], ab_id_pair[1], f'Antibody {ab_id_pair[0]} only contains background expression.'))
                    logging.warning(f"Skipping {ab_id_pair[0]}. Refer to antibody_errors.txt")
                    continue
                
                # If antibody does not exist, AND it is not completely background
                # Else, we want to do a lookup
                each_hit_results = _sci_uni(ab_id_pair)

            try:
                
                # Make request to Antibody Registry for additional info
                results_df = extra_ab_info(ab_id_pair[1])
                ab_name = str(results_df.loc[0]['ab_name'])
                ab_target = str(results_df.loc[0]['ab_target'])
                citation = str(results_df.loc[0]['proper_citation'])
                comments = str(results_df.loc[0]['comments'])
                vendor = str(results_df.loc[0]['vendor'])
                catalogNum = str(results_df.loc[0]['catalog_num']) # this will be a str varchar in db

                for hit_results in each_hit_results:
                    alias, clonalBool, host, uniprotID, cloneID, otherAliases = hit_results

                    cursor.callproc('insert_antibody',
                                    args=(ab_id_pair[1], alias, uniprotID, clonalBool, host, cloneID,
                                          ab_name, ab_target, citation, comments, vendor, catalogNum))
                    # Commit query changes
                    conn.commit()

                    for otherAlias in otherAliases:
                        cursor.callproc('insert_alias', args=(otherAlias, uniprotID))
                        conn.commit()
            except:
                raise Exception(f"Error with inserting antibody {ab_id_pair[1]} into database")
        
    # Log errors to an external file
    if errors:
        with open('antibody_errors.txt', 'a') as f:
            for error in errors:
                time = str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
                f.write(f"[{time}] {error}\n")
    
    if conn is not None:
        print("Disconnecting from database...\n")
        cursor.close()
        conn.close()

def convert_idBTO_readable(idBTO: str) -> str:
    """
    Converts a tissue ontology ID into its common tissue name

    Parameters:
        idBTO (str): tissue ontology ID
    
    Returns:
        tissueType (str): tissue name from tissue ontology ID
    """

    idBTO_params = {
        'q': idBTO,
        'ontology': 'bto',
        'local': 'true',
        'exact': 'true',
        'rows': 1,
        'start': 0
    }
    
    res = requests.get(EBI_BASE, params=idBTO_params)
    res_JSON = res.json()
    tissueType = res_JSON['response']['docs'][0]['label']
    
    return tissueType

def _connect_db_experiment(specs_csv: str, IPD) -> int:
    """
    Connects to a database and inserts information about an experiment

    Parameters:
        specs_csv (str): name of csv file containing a spreadsheet with
            information about the experiment and antibodies
        IPD (ImmunoPhenoData object): object containing protein data, 
            gene data, cell labels, and cell certainties
    
    Returns:
        id_result (int): the id of the experiment that was recently inserted
            into the database. Will be used in other procedures.
    """
    protein_counts = IPD._temp_protein

    # CHECK: Antibodies in specs spreadsheet must be in the protein data
    ab_id_pairs = _read_antibodies(specs_csv)
    ab_names = set([ab_id[0] for ab_id in ab_id_pairs])

    # differences = set(ab_names).difference(set(protein_counts.columns))
    # # If there are antibodies in the spreadsheet that are not in the data, raise exception
    # if len(differences) > 0:
    #     for ab in differences:
    #         logging.warning(f"ERROR: Antibody {ab} in spreadsheet not found in dataset")
        
    #     raise Exception("Mismatched antibody names between spreadsheet and protein data.")

    experiment = _read_experiment(specs_csv)

    expName = experiment['exp']
    expType = experiment['type']

    if 'cite' in expType.lower():
        expTypeEntry = 'CITE'
    elif 'reap' in expType.lower():
        expTypeEntry = 'REAP'
    elif 'ab' in expType.lower():
        expTypeEntry = 'AB'

    expPMID = experiment['pmid']
    expDOI = experiment['doi']
    expTissue = experiment['tissue']

    # Connect to database
    params = _config()
    print('Connecting to the MySQL database...')
    conn = mysql.connector.connect(**params)
    print("Connected to db\n")
    cursor = conn.cursor()

    try:
        print("Inserting experiment name:", expName)
        print("Inserting experiment type:", expTypeEntry)
        print("Inserting experiment PMID:", expPMID)
        print("Inserting experiment DOI:", expDOI)
        print("Inserting experiment tissue:", expTissue)

        # Make request to OLS database for tissue name
        tissue_name = convert_idBTO_readable(expTissue)

        check_exp_exists_query = """SELECT COUNT(*)
                                    FROM experiments
                                    WHERE experiments.nameExp=(%s)"""
        # Parameter must be converted from str to tuple
        cursor.execute(check_exp_exists_query, (expName, ))
        exp_exists_result = cursor.fetchone()[0]

        if exp_exists_result != 1:
            # If experiment doesn't exist, we can add it
            cursor.callproc("insert_experiment", 
                            args=(expName, expTypeEntry, expPMID, expDOI, expTissue, tissue_name))
            conn.commit()
        else:
            # If it does exist, do nothing
            pass

        sql_search_query = """SELECT experiments.idExperiment
                            FROM experiments
                            ORDER BY idExperiment DESC
                            LIMIT 0, 1"""

        cursor.execute(sql_search_query)
        id_result = cursor.fetchone()[0]
        print(f"Experiment ID: {id_result}")
    except:
        raise Exception("Error with inserting experiment")

    if conn is not None:
        print("Disconnecting from database...\n")
        cursor.close()
        conn.close()
    
    return id_result

def _experiment_idCell_map(idExperiment: int) -> dict:
    """
    Creates a dictionary using the database 'cells' table
    
    Parameters:
        idExperiment (int): experiment ID to retrieve all cells
            for that particular experiment
        
    Returns:
        cells_dict (dict): dictionary containing
            key: barcode tags for individual cells
            values: idCell, an auto-increment value for a given cell 
    """
    
    # Connect to database
    params = _config()
    print('Connecting to the MySQL database...')
    conn = mysql.connector.connect(**params)
    print("Connected to db\n")
    cursor = conn.cursor()
    
    cells_table_query = """SELECT cells.idCellOriginal, cells.idCell
                        FROM cells
                        WHERE cells.idExperiment = (%s);""" % (idExperiment)
    
    cells_df = pd.read_sql(sql=cells_table_query, con=conn)
    
    # Create mapping dict using columns "idCellOriginal" and "idCell"
    cells_dict = dict(zip(cells_df.idCellOriginal, cells_df.idCell))
    
    return cells_dict

def _connect_db_antigen_expression(idExperiment: int,
                                   specs_csv: str,
                                   IPD):
    """
    Connects to a database and populates the antigen_expression table

    Parameters:
        idExperiment (int): the id of the experiment inserted in the database
        specs_csv (str): name of csv file containing a spreadsheet with
            information about the experiment and antibodies
        IPD (ImmunoPhenoData object): object containing protein data, 
            gene data, cell labels, and cell certainties  
    """
    # Connect to database
    params = _config()
    print('Connecting to the MySQL database...')
    conn = mysql.connector.connect(**params)
    print("Connected to db\n")
    cursor = conn.cursor()
    cursor_idCell = conn.cursor()

    errors = []
    ab_lookup = _ab_dict(specs_csv)

    raw_counts = IPD.protein
    classified_counts = IPD.classified_filt
    normalized_counts = IPD.normalized_counts

    num_antibodies = len(normalized_counts.columns)

    print("Inserting into antigen-expression table...")

    # Iterate over columns, with antibody name and count values
    with _tqdm_output(tqdm(normalized_counts.items(), total=num_antibodies)) as tqdm_normalized:
        # Create a map for all cells for this given experiment
        # This will contain "idCellOriginal" : "idCell" as a dictionary entry pair
        cells_map = _experiment_idCell_map(idExperiment)

        for ab, counts in tqdm_normalized:
            # Get the corresponding antibody id for this antibody
            ab_id = ab_lookup[ab]

            # First check if this antibody is completely background       
            all_bg = all(x == 0 for x in (IPD.classified_filt.loc[:, ab]))
            if all_bg:
                # Skip this antibody
                errors.append((ab, ab_id, 'Antibody contains only background expression.'))
                logging.warning(f"Skipping {ab}. Refer to antigen_errors.txt")
                continue

            # CHECK: Does this antibody id exist in the database?
            check_ab_exists_query = """SELECT COUNT(*)
                                    FROM antibodies 
                                    WHERE antibodies.idAntibody=(%s)"""
            # Parameter must be converted from str to tuple
            cursor.execute(check_ab_exists_query, (ab_id, ))
            ab_exists_result = cursor.fetchone()[0]
            
            # If antibody exists, the value will be 1
            if ab_exists_result == 1:
                # Start parsing through each cell for this antibody
                for cell_name, value in counts.items():
                    # Get the matching idCell for this cell 
                    # idCell_query = """SELECT idCell
                    #                 FROM cells
                    #                 WHERE cells.idCellOriginal=(%s)
                    #                 AND cells.idExperiment=(%s)"""

                    # cursor_idCell.execute(idCell_query, (cell_name, idExperiment))
                    
                    # idCell_res = cursor_idCell.fetchone()
                    idCell = cells_map[cell_name]

                    # # Check if idCell exists in database
                    # if idCell_res is not None:
                    #     idCell = idCell_res[0]
                    # # If it doesn't, then skip this cell.
                    # else:
                    #     continue

                    # Grab the raw value of this cell-antibody match
                    ab_raw_val = int(raw_counts.loc[cell_name][ab])

                    # Grab the normalized value
                    ab_norm_val = float(normalized_counts.loc[cell_name][ab])

                    # Get the classification status 
                    ab_classification = int(classified_counts.loc[cell_name][ab])

                    # In the classified_csv, 0 means background cell, 1 is signal
                    # BUT in MySQL, 0 represents false (not background cell)
                    # We need to swap here
                    if ab_classification == 0:
                        background = True
                    elif ab_classification == 1:
                        background = False
                    
                    try:
                        cursor.callproc("insert_antigen_expression", 
                                        args=(idCell,
                                            ab_id,
                                            idExperiment,
                                            ab_raw_val,
                                            ab_norm_val,
                                            background))
                        conn.commit()
                    except:
                        raise Exception("Error with inserting antigen-expression value")

            # If it does not exist, it will be 0. Store in error log
            elif ab_exists_result == 0:
                errors.append((ab, ab_id, 'Antibody does not exist in database.'))
                logging.warning(f"Skipping {ab}. Refer to antigen_errors.txt")
                continue
    
    if conn is not None:
        print("Disconnecting from database...\n")
        cursor.close()
        cursor_idCell.close()
        conn.close()
    
    # Log errors to an external file
    if errors:
        with open('antigen_errors.txt', 'a') as f:
            f.write("\n")
            for error in errors:
                time = str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
                f.write(f"[{time}] {error}\n")

def convert_idCL_readable(idCL:str):
    """
    Converts a cell ontology ID into its common cell type name

    Parameters:
        idCL (str): cell ontology ID
    
    Returns:
        cellType (str): cell type name from cell ontology ID
    """

    idCL_params = {
        'q': idCL,
        'ontology': 'cl',
        'local': 'true',
        'exact': 'true',
        'rows': 1,
        'start': 0
    }
    
    res = requests.get(EBI_BASE, params=idCL_params)
    res_JSON = res.json()
    cellType = res_JSON['response']['docs'][0]
    
    return cellType

def _connect_db_cells(idExperiment: int,
                     specs_csv: str,
                     IPD):
    """
    Connects to a database and inserts information into the cells table
    and updates the antigen_expression table.

    Parameters:
        idExperiment (int): the id of the experiment inserted in the database
        specs_csv (str): name of csv file containing a spreadsheet with
            information about the experiment and antibodies
        IPD (ImmunoPhenoData object): object containing protein data, 
            gene data, cell labels, and cell certainties
    """
    # Connect to database
    params = _config()
    print('Connecting to the MySQL database...')
    conn = mysql.connector.connect(**params)
    print("Connected to db\n")
    cursor = conn.cursor()

    # Check if idExperiment is valid
    check_experiment_exists_query = """SELECT COUNT(*)
                                    FROM experiments
                                    WHERE experiments.idExperiment=(%s)"""
    cursor.execute(check_experiment_exists_query, (idExperiment, ))
    experiment_exists_result = cursor.fetchone()[0]
    if experiment_exists_result == 0:
        raise Exception("Error. idExperiment not found in database")

    # Check if IPD contains cell labels and certainties
    if IPD.raw_cell_labels is None or IPD.label_certainties is None:
        raise Exception("Error. No cell labels found. Run annotate_cells() first.")

    labels = IPD.raw_cell_labels["celltype"]
    deltas = IPD.label_certainties

    # Apply any filtering by delta or certainity value
    combined = pd.concat([labels, deltas], axis=1)

    # Apply same filtering to cells that were filtered out during normalization
    combined_norm = combined.loc[IPD.normalized_counts.index]

    print("Inserting new cells...")  
    # Insert cells that are in the combined singleR dataframe
    with _tqdm_output(tqdm(combined_norm.iterrows(), total=len(combined_norm))) as tqdm_cells:
        for index, row in tqdm_cells:
            check_cell_exists_query = """SELECT COUNT(*)
                                    FROM cells
                                    WHERE cells.idCellOriginal=(%s)
                                    AND cells.idExperiment=(%s)"""
            # Parameter must be converted from str to tuple
            cursor.execute(check_cell_exists_query, (index, idExperiment, ))
            cell_exists_result = cursor.fetchone()[0]
            if cell_exists_result == 0:
                try:
                    # Get readable idCL name and link
                    # idCL_info = convert_idCL_readable(row['labels'])
                    # cell_type_name = idCL_info['label']

                    # EDIT: We already have the readable name from annotate_cells, no need to do it here
                    cell_type_name = row['celltype']

                    cursor.callproc("insert_cell", args=(row['labels'], index, idExperiment, row['delta.next'],
                                                         cell_type_name))
                    conn.commit()
                except:
                    raise Exception("Error with inserting cell")
            else:
                continue
    
    # Remove cells from the database that are no longer in the combined singleR dataframe
    all_cells_query = """SELECT cells.idCellOriginal
                        FROM cells
                        WHERE cells.idExperiment=(%s);"""
    cursor.execute(all_cells_query, (idExperiment ,))
    db_cells = cursor.fetchall()
    db_cells_list = [cell[0] for cell in db_cells]

    # Cells in dataframe
    sr_cells = combined_norm.index

    # Find difference between two using XOR operator
    diff = set(db_cells_list) ^ set(sr_cells)

    # Remove these cells from the database
    # this will remove rows from the 'cells' table and 'antigen_expression' table
    print("Removing previous cells... ")
    with _tqdm_output(tqdm(diff, total=len(diff))) as tqdm_cells:
        remove_foreign_key_check = """SET FOREIGN_KEY_CHECKS=0"""
        cursor.execute(remove_foreign_key_check)
        conn.commit()

        for cell_name in tqdm_cells:
            remove_cell_query = """DELETE FROM cells
                                WHERE cells.idCellOriginal=(%s)
                                AND cells.idExperiment=(%s)
                                """
            try:
                # Parameter must be converted from str to tuple
                cursor.execute(remove_cell_query, (cell_name, idExperiment, ))
                conn.commit()
            except:
                raise Exception("Error removing past cells.")
        
        add_foreign_key_check = """SET FOREIGN_KEY_CHECKS=1"""
        cursor.execute(add_foreign_key_check)
        conn.commit()

    if conn is not None:
        print("Disconnecting from database...\n")
        cursor.close()
        conn.close()
    
    # Call connect_db_antigen_expression to fill antigen_expression table
    _connect_db_antigen_expression(idExperiment,
                                  specs_csv,
                                  IPD)

def link_ab_uniprot(ab_id_pair: list,
                    uniprotID: str):
    """
    Insert an antibody into the database using an explicit UniProtID

    Parameters:
        ab_id_pair (list[str, str]): antibody name and id. Ex: ['CD90', 'AB_123']
        uniprotID (str): accession ID from UniProtKB
    """
    # Connect to database
    params = _config()
    print('Connecting to the MySQL database...')
    conn = mysql.connector.connect(**params)
    print("Connected to db\n")
    cursor = conn.cursor()

    # Retrieve results from SciCrunch and UniProt
    hits = _sci_uni(ab_id_pair, uniprotID)
    try:
        for hit in hits:
            alias, clonalBool, host, uniprotID, cloneID, otherAliases = hit

            # Make request to Antibody Registry for additional info
            results_df = extra_ab_info(ab_id_pair[1])
            ab_name = str(results_df.loc[0]['ab_name'])
            ab_target = str(results_df.loc[0]['ab_target'])
            citation = str(results_df.loc[0]['proper_citation'])
            comments = str(results_df.loc[0]['comments'])
            vendor = str(results_df.loc[0]['vendor'])
            catalogNum = str(results_df.loc[0]['catalog_num']) # this will be a str varchar in db

            print(f"Manually inserting {ab_id_pair[0]}...")
            cursor.callproc('insert_antibody', 
                            args=(ab_id_pair[1], alias, uniprotID, clonalBool, host, cloneID,
                                  ab_name, ab_target, citation, comments, vendor, catalogNum))
            
            # Commit query changes
            conn.commit()

            for otherAlias in otherAliases:
                cursor.callproc('insert_alias', args=(otherAlias, uniprotID))
                conn.commit()
    except:
        raise Exception(f"Error manually inserting antibody {ab_id_pair[1]} into database")

    if conn is not None:
        print("Disconnecting from database...\n")
        cursor.close()
        conn.close()

def link_antigen(ab_id_pair: list,
                 idExperiment: int,
                 IPD):
    """
    Populate antigen_expression table for a single antibody

    Parameters:
        ab_id_pair (list[str, str]): antibody name and id. Ex: ['CD90', 'AB_123']
        idExperiment (int): the id of the experiment inserted in the database
        IPD (ImmunoPhenoData object): object containing protein data, 
            gene data, cell labels, and cell certainties  
    """
    # Connect to database
    params = _config()
    print('Connecting to the MySQL database...')
    conn = mysql.connector.connect(**params)
    print("Connected to db\n")
    cursor = conn.cursor()

    # Assumes link_ab_uniprot has been called already to fill in the missing antibody
    ab_name = ab_id_pair[0]
    ab_id = ab_id_pair[1]

    raw_counts = IPD.protein
    classified_counts = IPD.classified_filt
    normalized_counts = IPD.normalized_counts
    num_cells = len(normalized_counts.index)

    # CHECK: Does this antibody id exist in the database?
    check_ab_exists_query = """SELECT COUNT(*)
                            FROM antibodies 
                            WHERE antibodies.idAntibody=(%s)"""
    # Parameter must be converted from str to tuple
    cursor.execute(check_ab_exists_query, (ab_id, ))
    
    ab_exists_result = cursor.fetchone()[0]

    if ab_exists_result == 1:
        print(f"Manually inserting {ab_name} into antigen expression table...")
        ab_norm_counts = normalized_counts.loc[:, ab_name].items()
        # Look into normalized counts
        with _tqdm_output(tqdm(ab_norm_counts, total=num_cells)) as tqdm_normalized:
            # Create a map for all cells for this given experiment
            # This will contain "idCellOriginal" : "idCell" as a dictionary entry pair
            cells_map = _experiment_idCell_map(idExperiment)

            for cell_name, value in tqdm_normalized:
                # Get the matching idCell for this cell 
                # idCell_query = """SELECT idCell
                #                 FROM cells
                #                 WHERE cells.idCellOriginal=(%s)
                #                 AND cells.idExperiment=(%s)"""

                # cursor.execute(idCell_query, (cell_name, idExperiment))
                
                # idCell_res = cursor.fetchone()
                idCell = cells_map[cell_name]

                # # Check if idCell exists in database
                # if idCell_res is not None:
                #     idCell = idCell_res[0]
                # # If it doesn't, then skip this cell.
                # else:
                #     continue

                # Grab the raw value of this cell-antibody match
                ab_raw_val = int(raw_counts.loc[cell_name][ab_name])

                # Grab the normalized value
                ab_norm_val = float(value)

                # Get the classification status 
                ab_classification = int(classified_counts.loc[cell_name][ab_name])

                if ab_classification == 0:
                    background = True
                elif ab_classification == 1:
                    background = False
                
                try:
                    cursor.callproc("insert_antigen_expression", 
                                    args=(idCell,
                                        ab_id,
                                        idExperiment,
                                        ab_raw_val,
                                        ab_norm_val,
                                        background))
                    conn.commit()
                except:
                    raise Exception("Error with manually inserting antigen-expression value")

    elif ab_exists_result == 0:
        raise Exception(f"Unable to locate antibody ID:{ab_id} in database. Re-check parameters and call link_ab_uniprot()")

    if conn is not None:
        print("Disconnecting from database...\n")
        cursor.close()
        conn.close()

def load_csv_database(IPD):
    """
    Wrapper function to populate the entire database using a spreadsheet
    containing an experiment and antibodies along with an ImmunoPhenoData object
    from the ImmunoPheno package

    Parameters:
        specs_csv (str): name of csv file containing a spreadsheet with
            information about the experiment and antibodies
        IPD (ImmunoPhenoData object): object containing protein data, 
            gene data, cell labels, and cell certainties 
    """
    
    # Create tables and procedures
    _connect_db_tables()
    _connect_db_procedures()

    # Retrieve idExperiment after putting into database
    idExperiment = _connect_db_experiment(IPD._spreadsheet, IPD)
    _connect_db_antibody(IPD._spreadsheet, IPD)
    _connect_db_cells(idExperiment, IPD._spreadsheet, IPD)

def experiments() -> pd.DataFrame:
    """
    Returns all experiments in the current database
    """
    
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)
    
    experiments_query = """SELECT * FROM experiments;"""
    experiments_df = pd.read_sql(sql=experiments_query, con=conn)
    
    return experiments_df

def delete_experiment(idExperiment: int):
    params = _config()
    print('Connecting to the MySQL database...')
    conn = mysql.connector.connect(**params)
    print("Connected to db\n")
    cursor = conn.cursor()
    
    # Before starting, check if experiment even exists in database
    exp_exists_query = """SELECT COUNT(*)
                        FROM experiments
                        WHERE experiments.idExperiment = (%s)"""
    cursor.execute(exp_exists_query, (idExperiment, ))
    exists_result = cursor.fetchone()[0]
    
    if exists_result == 0:
        raise Exception(f"Error. Experiment {idExperiment} does not exist in the database.")
        
    # Remove foreign key check
    remove_foreign_key_check = """SET FOREIGN_KEY_CHECKS=0"""
    cursor.execute(remove_foreign_key_check)
    conn.commit()
    
    idExp_antibodies_query = """SELECT DISTINCT antigen_expression.idAntibody
                            FROM antigen_expression
                            WHERE antigen_expression.idExperiment = (%s)"""
    
    cursor.execute(idExp_antibodies_query, (idExperiment, ))
    all_abs = cursor.fetchall()
    all_abs_list = [ab[0] for ab in all_abs]
    
    ab_to_delete = []
    
    # Find which antibodies need to be removed from the 'antibodies' table
    for ab in all_abs_list:
        num_idExp_query = """SELECT COUNT(DISTINCT (antigen_expression.idExperiment))
                        FROM antigen_expression
                        WHERE antigen_expression.idAntibody = (%s)"""
        cursor.execute(num_idExp_query, (ab, ))
        num_idExp = cursor.fetchone()[0]
        
        if num_idExp == 1:
            ab_to_delete.append(ab)
    
    # Delete all rows from antigen_expression for this experiment (regardless of above)
    delete_antigen_expression_query = """DELETE FROM antigen_expression
                                    WHERE antigen_expression.idExperiment = (%s) """
    print(f"Deleting rows in 'antigen_expression' for experiment {idExperiment}...")
    cursor.execute(delete_antigen_expression_query, (idExperiment, ))
    conn.commit()
            
    # Delete all related information about the antibody from the 'antibodies' table
    for ab in ab_to_delete:
        idUniprot_query = """SELECT antigens.idUniProtKB
                        FROM antigens
                        WHERE antigens.idAntibody = (%s)"""
        cursor.execute(idUniprot_query, (ab, ))
        idUniprot = cursor.fetchone()[0]
    
        delete_aliases_query = """DELETE FROM aliases
                            WHERE aliases.idUniProtKB = (%s)"""
        print(f"Deleting rows in 'aliases' for antibody: {ab}...")
        cursor.execute(delete_aliases_query, (idUniprot, ))
        conn.commit()
        
        delete_antigens_query = """DELETE FROM antigens
                                WHERE antigens.idAntibody = (%s)"""
        print(f"Deleting rows in 'antigens' for antibody: {ab}...")
        cursor.execute(delete_antigens_query, (ab, ))
        conn.commit()
        
        delete_proteins_query = """DELETE FROM proteins
                                WHERE proteins.idUniProtKB = (%s) """
        print(f"Deleting rows in 'proteins' for antibody: {ab}...")
        cursor.execute(delete_proteins_query, (idUniprot, ))
        conn.commit()
        
        delete_antibodies_query = """DELETE FROM antibodies
                                WHERE antibodies.idAntibody = (%s)"""
        print(f"Deleting rows in 'antibodies' for antibody: {ab}...")
        cursor.execute(delete_antibodies_query, (ab, ))
        conn.commit()
    
    # Find all cell types for this experiment
    idExp_celltypes_query = """SELECT DISTINCT cells.idCL
                            FROM cells
                            WHERE cells.idExperiment = (%s) """
    
    cursor.execute(idExp_celltypes_query, (idExperiment, ))
    all_idCLs = cursor.fetchall()
    all_idCLs_list = [idCL[0] for idCL in all_idCLs]
    
    idCL_to_delete = []
    
    # Find which idCLs need to be removed
    for idCL in all_idCLs_list:
        num_idExp_ct_query = """SELECT COUNT(DISTINCT (cells.idExperiment))
                            FROM cells
                            WHERE cells.idCL = (%s)"""
        cursor.execute(num_idExp_ct_query, (idCL, ))
        num_idExp = cursor.fetchone()[0]
        
        if num_idExp == 1:
            idCL_to_delete.append(idCL)
    
    # Before deleting cell types, delete all cells for this experiment
    delete_cells_query = """DELETE FROM cells
                        WHERE cells.idExperiment = (%s)"""
    print(f"Deleting rows in 'cells' for experiment {idExperiment}...")
    cursor.execute(delete_cells_query, (idExperiment, ))
    conn.commit()
    
    # Delete cell types 
    for idCL in idCL_to_delete:
        delete_idCL_query = """DELETE FROM cell_types
                            WHERE cell_types.idCL = (%s)"""
        print(f"Deleting rows in 'cell_types' for cell type {idCL}...")
        cursor.execute(delete_idCL_query, (idCL, ))
        conn.commit()
            
    # Find the current tissue of the experiment
    idExp_tissue_query = """SELECT experiments.idBTO
                        FROM experiments
                        WHERE experiments.idExperiment = (%s);"""
    cursor.execute(idExp_tissue_query, (idExperiment, ))
    idExp_tissue = cursor.fetchone()[0]
    
    # Find the number of experiments that use this tissue
    num_tissue_idExp_query = """SELECT COUNT(*)
                                FROM experiments
                                WHERE experiments.idBTO = (%s)"""
    cursor.execute(num_tissue_idExp_query, (idExp_tissue, ))
    num_tissue = cursor.fetchone()[0]
    
    # Remove the experiment from the experiment table
    delete_experiment_query = """DELETE FROM experiments
                            WHERE experiments.idExperiment = (%s)"""
    cursor.execute(delete_experiment_query, (idExperiment, ))
    conn.commit()
    
    # If there was only one experiment, remove the tissue
    if num_tissue == 1:
        delete_tissue_query = """DELETE FROM tissues
                            WHERE tissues.idBTO = (%s)"""
        print(f"Deleting rows in 'tissues' for tissue {idExp_tissue}...")
        cursor.execute(delete_tissue_query, (idExp_tissue, ))
        conn.commit()
        
    print(f"Deleting row in 'experiments' for experiment {idExperiment}...") 
    
    # Re-add foreign key check
    add_foreign_key_check = """SET FOREIGN_KEY_CHECKS=1"""
    cursor.execute(add_foreign_key_check)
    conn.commit()
        
    if conn is not None:
        print("Disconnecting from database...\n")
        cursor.close()
        conn.close()