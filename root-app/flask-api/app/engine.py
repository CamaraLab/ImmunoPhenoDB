import warnings
import logging

import configparser
import mysql.connector # mysql-connector-python

import requests

import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

SCI_CRUNCH_BASE = "http://www.scicrunch.org"
SCI_RRID_ENDPOINT = "/resolver/"
SCI_FILE = ".json"

EBI_BASE = "https://www.ebi.ac.uk/ols/api/search"

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

########################### Functions for get_antibodies ########################
def convert_ab_readable(ab_id:str):
    res = requests.get(SCI_CRUNCH_BASE + SCI_RRID_ENDPOINT + ab_id + SCI_FILE)
    res_JSON = res.json()
    ab_name = (res_JSON['hits']['hits'][0]
                    ['_source']['antibodies']['primary'][0]
                    ['targets'][0]['name'])
    return ab_name

def get_unique_idCLs():
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()
    
    idCLs_query = """SELECT unique_idCLs.idCL
                     FROM unique_idCLs;"""
    
    cursor.execute(idCLs_query)
    res = cursor.fetchall()
    idCLs_list = [idCL[0] for idCL in res]
    
    if conn is not None:
        cursor.close()
        conn.close()
    
    return idCLs_list

def get_antibodies_idCL(id_CLs: list,
                        idBTO: list = None,
                        idExperiment: list = None) -> list:
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()
    
    # Dynamically generate placeholders for parameterized query
    idCL_placeholders = ','.join(['%s'] * len(id_CLs))
    
    # No additional filtering
    if idBTO is None and idExperiment is None:
        ab_query = """SELECT DISTINCT antigen_expression.idAntibody
                    FROM antigen_expression
                    WHERE antigen_expression.idCell IN 
                        (SELECT cells.idCell
                        FROM cells
                        WHERE cells.idCL IN (%s) );""" % (idCL_placeholders)
        cursor.execute(ab_query, (id_CLs))
    
    # Filtering based on tissue type (idBTO)
    elif idBTO is not None and idExperiment is None:
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
        ab_idBTO_query = """SELECT DISTINCT antigen_expression.idAntibody
                        FROM antigen_expression
                        WHERE antigen_expression.idCell IN 
                            (SELECT cells.idCell
                            FROM cells
                            WHERE cells.idCL IN (%s) AND
                            cells.idExperiment IN 
                                (SELECT experiments.idExperiment
                                FROM experiments
                                WHERE experiments.idBTO IN (%s) ));""" % (idCL_placeholders, idBTO_placeholders)
        temp_idCLs = id_CLs.copy()
        temp_idCLs.extend(idBTO)

        cursor.execute(ab_idBTO_query, (temp_idCLs))
    
    # Filtering based on experiment id (idExperiment)
    elif idBTO is None and idExperiment is not None:
        idExperiment_placeholders = ','.join(['%s'] * len(idExperiment))
        ab_idExperiment_query = """SELECT DISTINCT antigen_expression.idAntibody
                                FROM antigen_expression
                                WHERE antigen_expression.idCell IN 
                                    (SELECT cells.idCell
                                    FROM cells
                                    WHERE cells.idCL IN (%s) AND
                                    cells.idExperiment IN (%s) );""" % (idCL_placeholders, idExperiment_placeholders)
        temp_idCLs = id_CLs.copy()
        temp_idCLs.extend(idExperiment)
        cursor.execute(ab_idExperiment_query, (temp_idCLs))
    
    # Filtering based on both idBTO and idExperiment
    elif idBTO is not None and idExperiment is not None:
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
        idExperiment_placeholders = ','.join(['%s'] * len(idExperiment))
        print('do we enter?')
        ab_idBTO_idExp_query = """SELECT DISTINCT antigen_expression.idAntibody
                                FROM antigen_expression
                                WHERE antigen_expression.idCell IN 
                                (SELECT cells.idCell
                                FROM cells
                                WHERE cells.idCL IN (%s) AND
                                    cells.idExperiment IN (%s) AND
                                    cells.idExperiment IN 
                                (SELECT experiments.idExperiment
                                FROM experiments
                                WHERE experiments.idBTO IN (%s)));""" % (idCL_placeholders, idExperiment_placeholders, idBTO_placeholders)
        temp_idCLs = id_CLs.copy()
        temp_idCLs.extend(idExperiment)
        temp_idCLs.extend(idBTO)
        print(temp_idCLs)
        cursor.execute(ab_idBTO_idExp_query, (temp_idCLs))
    
    res = cursor.fetchall()
    ab_list = [ab[0] for ab in res]
    
    if conn is not None:
        cursor.close()
        conn.close()
    
    return ab_list

def get_cells_ab_idCL(ab_id: str,
                      id_CLs: list):
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)

    # Dynamically generate placeholders for parameterized query
    idCL_placeholders = ','.join(['%s'] * len(id_CLs))
    
    with_query = """SELECT cells.idCell, cells.idCellOriginal, antigen_expression.normValue, cells.idExperiment
                FROM antigen_expression
                INNER JOIN cells ON antigen_expression.idCell = cells.idCell
                WHERE antigen_expression.idCell IN 
                    (SELECT cells.idCell
                    FROM cells
                    WHERE cells.idCell IN 
                        (SELECT antigen_expression.idCell
                        FROM antigen_expression
                        WHERE antigen_expression.idAntibody = (%s) ) AND
                        cells.idCL IN (%s) ) AND
                    antigen_expression.idAntibody = (%s);""" % ('%s', idCL_placeholders, '%s')
    
    without_query = """SELECT cells.idCell, cells.idCellOriginal, antigen_expression.normValue, cells.idExperiment
                    FROM antigen_expression
                    INNER JOIN cells ON antigen_expression.idCell = cells.idCell
                    WHERE antigen_expression.idCell IN 
                        (SELECT cells.idCell
                        FROM cells
                        WHERE cells.idCell IN 
                            (SELECT antigen_expression.idCell
                            FROM antigen_expression
                            WHERE antigen_expression.idAntibody = (%s) ) AND
                            cells.idCL NOT IN (%s) ) AND
                        antigen_expression.idAntibody = (%s);""" % ('%s', idCL_placeholders, '%s')
    
    parameters = [ab_id]
    parameters.extend(id_CLs)
    parameters.append(ab_id)
    
    cells_with_idCL_df = pd.read_sql(sql=with_query, params=parameters, con=conn)
    cells_with_idCL_df['cellType'] = 0
    
    cells_without_idCL_df = pd.read_sql(sql=without_query, params=parameters, con=conn)
    cells_without_idCL_df['cellType'] = 1

    # Combine both dataframes
    combined = pd.concat([cells_with_idCL_df, cells_without_idCL_df], ignore_index=True)
    
    if conn is not None:
        conn.close()
        
    return combined

def retrieve_lmm(summary_table):
    descr = summary_table.summary().tables[0]
    
    # Return top half as a dictionary
    descr_left = dict(zip(descr[0], descr[1]))
    descr_right = dict(zip(descr[2], descr[3]))
    descr_dict = {**descr_left, **descr_right}
    
    # Return bottom half as a dataframe
    core = summary_table.summary().tables[1]
        
    return (descr_dict, core)

def run_lmm(id_CLs: list,
            idBTO: list = None,
            idExperiment: list = None) -> dict:
    
    res_dict = {}
    
    # Retrieve unique antibodies from idCLs
    ab_db = get_antibodies_idCL(id_CLs, idBTO, idExperiment)
    
    # Perform Linear Mixed Models for each antibody
    for ab in ab_db:
        
        # Cells associated with each antibody
        comb_df = get_cells_ab_idCL(ab, id_CLs)
        
        # Run LMM
        lmm = smf.mixedlm("normValue ~ cellType", comb_df, groups=comb_df["idExperiment"])
        summary_table = lmm.fit(method=["bfgs"])
        
        # Retrieve results (tuple)
        lmm_res = retrieve_lmm(summary_table)
        
        # Add in p-values manually (with increased precision)
        p_values = summary_table.pvalues
        lmm_res[1].loc['Intercept']['P>|z|'] = p_values[0]
        lmm_res[1].loc['cellType']['P>|z|'] = p_values[1]
        
        res_dict[ab] = lmm_res
        
    return res_dict

def get_ab_bg_cells(ab_id, node_fam_dict):
    # Temporarily combine all keys and values for a total list of idCLs
    parents = list(node_fam_dict.keys())
    children = [item for sublist in list(node_fam_dict.values()) for item in sublist]
    family = parents + children
    
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)

    # Dynamically generate placeholders for parameterized query
    idCL_placeholders = ','.join(['%s'] * len(family))
    
    with_query = """SELECT antigen_expression.idAntibody, cells.idCL, COUNT(antigen_expression.background) as total_cells, 
                        SUM(IF(antigen_expression.background=1, 1, 0)) as bg_cells,
                        SUM(IF(antigen_expression.background=1, 1, 0))/COUNT(antigen_expression.background) as proportion
                    FROM antigen_expression
                    INNER JOIN cells ON antigen_expression.idCell = cells.idCell
                    WHERE antigen_expression.idCell IN (SELECT cells.idCell
                                                        FROM cells
                                                        WHERE cells.idCell IN (SELECT antigen_expression.idCell
                                                                                FROM antigen_expression
                                                                                WHERE antigen_expression.idAntibody = (%s) ) AND
                                                        cells.idCL IN (%s) ) AND
                    antigen_expression.idAntibody = (%s)
                    GROUP BY cells.idCL;""" % ('%s', idCL_placeholders, '%s')
    
    without_query = """SELECT COUNT(*) as total, SUM(IF(antigen_expression.background=1, 1, 0)) as bg_cells,
                        SUM(IF(antigen_expression.background=1, 1, 0))/COUNT(*) as proportion
                    FROM antigen_expression
                    INNER JOIN cells ON antigen_expression.idCell = cells.idCell
                    WHERE antigen_expression.idCell IN (SELECT cells.idCell
                                                        FROM cells
                                                        WHERE cells.idCell IN (SELECT antigen_expression.idCell
                                                                                FROM antigen_expression
                                                                                WHERE antigen_expression.idAntibody = (%s) ) AND
                                                        cells.idCL IN (%s) ) AND
                    antigen_expression.idAntibody = (%s);""" % ('%s', idCL_placeholders, '%s')
    
    parameters = [ab_id]
    parameters.extend(family)
    parameters.append(ab_id)

    prop_cells_with_idCL_df = pd.read_sql(sql=with_query, params=parameters, con=conn)
    
    # Set our index to be the idCL
    new_df = prop_cells_with_idCL_df.set_index('idCL')
    df_idCLs = list(prop_cells_with_idCL_df['idCL'])
    
    # If the parents have 0 cells? Well, it's set to 0 initially
    aggregate_df = pd.DataFrame(data=0, index=node_fam_dict.keys(), columns=['total_cells', 'bg_cells'])
    
    # We want to aggregate idCLs that are parent-descendant related
    for key, value in node_fam_dict.items():

        # If our parent key has a value in the df
        if key in df_idCLs:
            # Add their total cells 
            aggregate_df.loc[key]['total_cells'] += new_df.loc[key]['total_cells']
            # Add their background cells
            aggregate_df.loc[key]['bg_cells'] += new_df.loc[key]['bg_cells']
            
        for child in value:
            if child in df_idCLs:
                # Add their total cells to their parent's total_cells column
                aggregate_df.loc[key]['total_cells'] += new_df.loc[child]['total_cells']
                # Add their background cells
                aggregate_df.loc[key]['bg_cells'] += new_df.loc[child]['bg_cells']
                
        # WARNING: if there are 0 cells for a given idCL for an antibody
        if aggregate_df.loc[key]['total_cells'] == 0:
            logging.warning(f"No cells found for type {key} in {ab_id}. "
                           "Background percentage will be set to NaN. ")
    
    # Calculate the proportions of background cells into a new column
    aggregate_df['proportion'] = aggregate_df['bg_cells']/aggregate_df['total_cells']

    # Separate dataframe for cells not belonging to idCLs
    prop_cells_without_idCL_df = pd.read_sql(sql=without_query, params=parameters, con=conn)

    return (aggregate_df, prop_cells_without_idCL_df)

def ab_lmm_table(lmm_results: dict, node_fam_dict: dict):
    # Find antibodies for index
    antibodies = list(lmm_results.keys())
    
    cellType_coeff = []
    stdError = []
    p_values = []
    
    # Find all parent idCLs
    parents = list(node_fam_dict.keys())
    
    # Initialize dictionary with empty lists for each parent
    idCL_dict = {k: [] for k in parents}
    # Last row will be for all cells not with our idCLs
    idCL_dict['not_idCLs'] = []
    
    for ab, value in lmm_results.items():
        df = value[1]
        cellType_coeff.append(df.loc['cellType']['Coef.'])
        stdError.append(df.loc['cellType']['Std.Err.'])
        p_values.append(df.loc['cellType']['P>|z|'])
        
        aggregate_df, without_idCLs_df = get_ab_bg_cells(ab, node_fam_dict)
        
        # For each parent idCL, add the proportions
        for key in parents:
            value = aggregate_df.loc[key]['proportion']

            # Find appropriate key
            idCL_dict[key].append(value)
        
        # Add last column in
        idCL_dict['not_idCLs'].append(without_idCLs_df['proportion'].values[0])

    # Apply Benjamini Hochberg correction
    q_values = multipletests(p_values, method='fdr_bh')[1]
    
    table_dict = {'coeff': cellType_coeff,
                'stderr': stdError, 
                'p_val': p_values,
                'q_val': q_values}
    
    entire_dict = {**table_dict, **idCL_dict}

    final_df = pd.DataFrame(index=antibodies, data=entire_dict)
    
    # Add readable column of antibody names
    ab_conv = []
    for ab in antibodies:
        ab_name = convert_ab_readable(ab)
        ab_conv.append(ab_name)
    
    # Add this column into df
    final_df.insert(0, 'target', ab_conv)

    return final_df

def get_Abs(node_fam_dict: dict,
            idBTO: list = None,
            idExperiment: list = None) -> pd.DataFrame:
    
    # Collapse dict into one list
    parents = list(node_fam_dict.keys())
    children = [item for sublist in list(node_fam_dict.values()) for item in sublist]
    family = parents + children

    # Run linear mixed model with given idCLs and their descendants
    lmm_results = run_lmm(family, idBTO, idExperiment)
    
    # Create a results table for each antibody using results from lmm
    results_table = ab_lmm_table(lmm_results, node_fam_dict)
    
    return results_table

def plot_antibodies(ab_ids: list, 
                    id_CLs: list):
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)
  
    # Dynamically generate placeholders for parameterized query
    ab_placeholders = ','.join(['%s'] * len(ab_ids))
    idCL_placeholders = ','.join(['%s'] * len(id_CLs)) # find descendants on jupyter first, then send them here
    
    with_query = """SELECT cells.idCell, cells.idCellOriginal, 
                antigen_expression.normValue, cells.idExperiment, antigen_expression.idAntibody, cells.idCL
                FROM antigen_expression
                INNER JOIN cells ON antigen_expression.idCell = cells.idCell
                WHERE antigen_expression.idCell IN 
                    (SELECT cells.idCell
                    FROM cells
                    WHERE cells.idCell IN 
                        (SELECT antigen_expression.idCell
                        FROM antigen_expression
                        WHERE antigen_expression.idAntibody IN (%s) ) AND
                        cells.idCL IN (%s) ) AND
                    antigen_expression.idAntibody IN (%s);""" % (ab_placeholders, idCL_placeholders, ab_placeholders)
    
    parameters = ab_ids.copy()
    parameters.extend(id_CLs)
    parameters.extend(ab_ids)
    
    cells_with_idCL_df = pd.read_sql(sql=with_query, params=parameters, con=conn)
    
    if conn is not None:
        conn.close()
    
    return cells_with_idCL_df

########################## Functions for get_celltype ########################
def convert_idCL_readable(idCL:str):
    idCL_params = {
        'q': idCL,
        'exact': 'true',
        'ontology': 'cl',
        'fieldList': 'label',
        'rows': 1,
        'start': 0
    }
    
    res = requests.get(EBI_BASE, params=idCL_params)
    res_JSON = res.json()
    cellType = res_JSON['response']['docs'][0]['label']
    
    return cellType

def get_idCLs_ab_ids(ab_id: str,
                     idBTO: list = None, 
                     idExperiment: list = None):
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()
    
    # No additional filtering
    if idBTO is None and idExperiment is None:
        idCL_query = """SELECT DISTINCT cells.idCL
                    FROM cells
                    WHERE cells.idCell IN (SELECT antigen_expression.idCell
                    FROM antigen_expression
                    WHERE antigen_expression.idAntibody = (%s) );""" % ('%s')
        
        cursor.execute(idCL_query, (ab_id, ))
    
    # Filtering based on idBTO
    elif idBTO is not None and idExperiment is None:
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))

        idCL_idBTO_query = """SELECT DISTINCT cells.idCL 
                            FROM cells 
                            WHERE cells.idExperiment IN 
                                (SELECT experiments.idExperiment 
                                FROM experiments 
                                WHERE experiments.idBTO IN (%s) ) AND
                            cells.idCell IN
                                (SELECT antigen_expression.idCell
                                FROM antigen_expression
                                WHERE antigen_expression.idAntibody IN (%s) );""" % (idBTO_placeholders, '%s')
        temp_ab = idBTO.copy()
        temp_ab.append(ab_id)
        
        cursor.execute(idCL_idBTO_query, (temp_ab))
    
    # Filtering based on idExperiment
    elif idBTO is None and idExperiment is not None:
        idExperiment_placeholders = ','.join(['%s'] * len(idExperiment))
        idCL_idExp_query = """SELECT DISTINCT cells.idCL 
                            FROM cells 
                            WHERE cells.idExperiment IN (%s) AND
                            cells.idCell IN 
                                (SELECT antigen_expression.idCell
                                FROM antigen_expression
                                WHERE antigen_expression.idAntibody IN (%s) );""" % (idExperiment_placeholders, '%s')
        temp_ab = idExperiment.copy()
        temp_ab.append(ab_id)
        
        cursor.execute(idCL_idExp_query, (temp_ab))
    
    # Filtering based on idBTO and idExperiment
    elif idBTO is not None and idExperiment is not None:
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
        idExperiment_placeholders = ','.join(['%s'] * len(idExperiment))
        idCL_idBTO_idExp_query = """SELECT DISTINCT cells.idCL 
                                FROM cells 
                                WHERE cells.idExperiment IN 
                                    (SELECT experiments.idExperiment 
                                    FROM experiments 
                                    WHERE experiments.idBTO IN (%s) ) AND
                                cells.idExperiment IN (%s) AND
                                cells.idCell IN 
                                    (SELECT antigen_expression.idCell
                                    FROM antigen_expression
                                    WHERE antigen_expression.idAntibody IN (%s));""" % (idBTO_placeholders, idExperiment_placeholders, '%s')
        temp_ab = idBTO.copy()
        temp_ab.extend(idExperiment)
        temp_ab.append(ab_id)
        cursor.execute(idCL_idBTO_idExp_query, (temp_ab))
        
    res = cursor.fetchall()
    
    idCL_list = [idCL[0] for idCL in res]
    
    if conn is not None:
        cursor.close()
        conn.close()
    
    return idCL_list

def get_cells_by_type(ab_id: str, id_CL: str):
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)

    with_idCL_query = """SELECT cells.idCell, cells.idCellOriginal, antigen_expression.normValue, cells.idExperiment, antigen_expression.idAntibody
                        FROM antigen_expression
                        INNER JOIN cells ON antigen_expression.idCell = cells.idCell
                        WHERE antigen_expression.idCell IN 
                        (SELECT cells.idCell
                            FROM cells
                            WHERE cells.idCell IN (SELECT antigen_expression.idCell
                                                    FROM antigen_expression
                                                    WHERE antigen_expression.idAntibody = (%s) ) AND
                            cells.idCL = (%s) AND
                        antigen_expression.idAntibody = (%s));""" % ('%s', '%s', '%s')
    
    without_idCL_query = """SELECT cells.idCell, cells.idCellOriginal, antigen_expression.normValue, cells.idExperiment, antigen_expression.idAntibody
                        FROM antigen_expression
                        INNER JOIN cells ON antigen_expression.idCell = cells.idCell
                        WHERE antigen_expression.idCell IN 
                        (SELECT cells.idCell
                            FROM cells
                            WHERE cells.idCell IN (SELECT antigen_expression.idCell
                                                    FROM antigen_expression
                                                    WHERE antigen_expression.idAntibody = (%s) ) AND
                            cells.idCL != (%s) AND
                        antigen_expression.idAntibody = (%s));""" % ('%s', '%s', '%s')
    
    parameters = [ab_id]
    parameters.append(id_CL)
    parameters.append(ab_id)
    
    cells_with_idCL_df = pd.read_sql(sql=with_idCL_query, params=parameters, con=conn)
    cells_with_idCL_df['cellType'] = 0
    
    cells_without_idCL_df = pd.read_sql(sql=without_idCL_query, params=parameters, con=conn)
    cells_without_idCL_df['cellType'] = 1

    # Combine both dataframes
    combined = pd.concat([cells_with_idCL_df, cells_without_idCL_df], ignore_index=True)
    # print(combined)
    
    if conn is not None:
        conn.close()
        
    return combined

def run_lmm_cell_type(ab_id: str,
                      idCL_db: list) -> dict:
    lmm_results_dict = {}
    
    # # Retrieve unique idCLs from antibody
    # idCL_db = get_idCLs_ab_ids(ab_id, idBTO, idExperiment)

    # Perform Linear Mixed Models for idCL
    for idCL in idCL_db:
        # Cells associated with each antibody
        comb_df = get_cells_by_type(ab_id, idCL)
        
        # Run LMM
        lmm = smf.mixedlm("normValue ~ cellType", comb_df, groups=comb_df["idExperiment"])
        summary_table = lmm.fit(method=["bfgs"])
        
        # Retrieve results (tuple)
        lmm_res = retrieve_lmm(summary_table)
        
        # Add in p-values manually (with increased precision)
        p_values = summary_table.pvalues
        lmm_res[1].loc['Intercept']['P>|z|'] = p_values[0]
        lmm_res[1].loc['cellType']['P>|z|'] = p_values[1]
        
        lmm_results_dict[idCL] = lmm_res
        
    return lmm_results_dict

def get_idCL_bg_cells(ab_id, id_CL):
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)

    with_query = """SELECT antigen_expression.idAntibody, cells.idCL, 
                        COUNT(antigen_expression.background) as total_cells, 
                        SUM(IF(antigen_expression.background=1, 1, 0)) as bg_cells,
                        SUM(IF(antigen_expression.background=1, 1, 0))/COUNT(antigen_expression.background) as proportion
                    FROM antigen_expression
                    INNER JOIN cells ON antigen_expression.idCell = cells.idCell
                    WHERE antigen_expression.idCell IN (SELECT cells.idCell
                                                        FROM cells
                                                        WHERE cells.idCell IN (SELECT antigen_expression.idCell
                                                                                FROM antigen_expression
                                                                                WHERE antigen_expression.idAntibody = (%s) ) AND
                                                        cells.idCL = (%s) AND
                        antigen_expression.idAntibody = (%s) )
                    GROUP BY antigen_expression.idAntibody, cells.idCL;""" % ('%s', '%s', '%s')
    
    parameters = [ab_id]
    parameters.append(id_CL)
    parameters.append(ab_id)

    prop_cells_with_idCL_df = pd.read_sql(sql=with_query, params=parameters, con=conn)

    # Set our index to be the idCL
    new_df = prop_cells_with_idCL_df.set_index('idAntibody')
    
    if new_df.loc[ab_id]['total_cells'] == 0:
        logging.warning(f"No cells found for type {id_CL} in {ab_id}. "
                       "Proportion will be set to NaN")
    
    new_df['proportion'] = new_df['bg_cells']/new_df['total_cells']
  
    return new_df

def ab_lmm_table_cell_type(lmm_results: dict, ab_id: list):
    # Find idCLs for index
    idCLs = list(lmm_results.keys())
    
    cellType_coeff = []
    stdError = []
    p_values = []
    idCL_dict = {}
    
    # Last row will be for background cell percentage with that idCL and antibody
    idCL_dict['background'] = []
   
    for idCL, value in lmm_results.items():
        df = value[1]
        cellType_coeff.append(df.loc['cellType']['Coef.'])
        stdError.append(df.loc['cellType']['Std.Err.'])
        p_values.append(df.loc['cellType']['P>|z|'])

        with_idCLs_df = get_idCL_bg_cells(ab_id, idCL)

        # Add last column in
        idCL_dict['background'].append(with_idCLs_df['proportion'].values[0])
    
    # Apply Benjamini Hochberg correction
    q_values = multipletests(p_values, method='fdr_bh')[1]
    
    table_dict = {'coeff': cellType_coeff,
                'stderr': stdError, 
                'p_val': p_values,
                'q_val': q_values}

    entire_dict = {**table_dict, **idCL_dict}

    final_df = pd.DataFrame(index=idCLs, data=entire_dict)
    
    # Add cellType column of names from id
    idCL_conv = []
    for idCL in idCLs:
        idCL_name = convert_idCL_readable(idCL)
        idCL_conv.append(idCL_name)
    
    # Add column to final df
    final_df.insert(0, 'cellType', idCL_conv)

    return final_df

def get_cellType(ab_ids:list,
                  idBTO: list = None, 
                  idExperiment: list = None):
    # Dictionary to store all LMM tables
    final_dict = {}
    
    for ab in ab_ids: 
        # Find all unique idCLs for this ab (this is done in run_lmm_cell_type)
        idCLs = get_idCLs_ab_ids(ab, idBTO, idExperiment)

        # If there are no idCLs, it means that no cells were found for this ab
        if len(idCLs) == 0:
            logging.warning(f"No cells found for {ab}. "
                             f"Skipping {ab}")
            continue
        
        idCL_lmm_dict = run_lmm_cell_type(ab, idCLs)

        res_df = ab_lmm_table_cell_type(idCL_lmm_dict, ab)

        # Add result to dict
        final_dict[ab] = res_df
    
    return final_dict

def plot_celltypes(ab_id: str, 
                   id_CLs: list):
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)
    
    idCL_placeholders = ','.join(['%s'] * len(id_CLs))

    with_idCL_query = """SELECT cells.idCell, cells.idCellOriginal, 
                                antigen_expression.normValue, cells.idExperiment, antigen_expression.idAntibody, cells.idCL
                        FROM antigen_expression
                        INNER JOIN cells ON antigen_expression.idCell = cells.idCell
                        WHERE antigen_expression.idCell IN 
                        (SELECT cells.idCell
                            FROM cells
                            WHERE cells.idCell IN (SELECT antigen_expression.idCell
                                                    FROM antigen_expression
                                                    WHERE antigen_expression.idAntibody = (%s) ) AND
                            cells.idCL IN (%s) AND
                        antigen_expression.idAntibody = (%s));""" % ('%s', idCL_placeholders, '%s')
    
    parameters = [ab_id]
    parameters.extend(id_CLs)
    parameters.append(ab_id)
    
    cells_with_idCL_df = pd.read_sql(sql=with_idCL_query, params=parameters, con=conn)
    
    if conn is not None:
        conn.close()
        
    return cells_with_idCL_df

########################## Functions for get_experiment ########################
def get_exp(ab_id: list, 
            idCL: list = None,
            idBTO: list = None): #function assumes we found the family idCLs from class
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()
    
    ab_placeholders = ','.join(['%s'] * len(ab_id))
    
    if idCL is None and idBTO is None:
        ab_query = """SELECT *
                    FROM experiments
                    WHERE experiments.idExperiment IN 
                        (SELECT DISTINCT antigen_expression.idExperiment
                         FROM antigen_expression
                         WHERE antigen_expression.idAntibody IN (%s));""" % (ab_placeholders)
        
        parameters = ab_id.copy()
        exp_df = pd.read_sql(sql=ab_query, params=parameters, con=conn)

    elif idCL is not None and idBTO is None:
        idCLs_placeholders = ','.join(['%s'] * len(idCL))
    
        idCL_query = """SELECT * 
                        FROM experiments
                        WHERE experiments.idExperiment IN 
                            (SELECT DISTINCT cells.idExperiment
                             FROM cells
                             WHERE cells.idCell IN (SELECT antigen_expression.idCell
                                                    FROM antigen_expression
                                                    WHERE antigen_expression.idAntibody IN (%s) ) 
                            AND cells.idCL IN (%s) );""" % (ab_placeholders, idCLs_placeholders)
        parameters = ab_id.copy()
        parameters.extend(idCL)
        exp_df = pd.read_sql(sql=idCL_query, params=parameters, con=conn)
        
    elif idCL is None and idBTO is not None:
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
        
        BTO_query = """SELECT * 
                        FROM experiments
                        WHERE experiments.idExperiment IN 
                            (SELECT DISTINCT cells.idExperiment
                                FROM cells
                                WHERE cells.idCell IN 
                                (SELECT antigen_expression.idCell
                                    FROM antigen_expression
                                    WHERE antigen_expression.idAntibody IN (%s) ) 
                        AND experiments.idBTO IN (%s) );""" % (ab_placeholders, idBTO_placeholders)
        parameters = ab_id.copy()
        parameters.extend(idBTO)
        exp_df = pd.read_sql(sql=BTO_query, params=parameters, con=conn)                 

    elif idCL is not None and idBTO is not None:
        idCLs_placeholders = ','.join(['%s'] * len(idCL))
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
    
        idCL_BTO_query = """SELECT * 
                            FROM experiments
                            WHERE experiments.idExperiment IN 
                                (SELECT DISTINCT cells.idExperiment
                                 FROM cells
                                 WHERE cells.idCell IN 
                                    (SELECT antigen_expression.idCell
                                     FROM antigen_expression
                                     WHERE antigen_expression.idAntibody IN (%s) ) 
                                AND cells.idCL IN (%s) 
                            AND experiments.idBTO IN (%s) );""" % (ab_placeholders, idCLs_placeholders, idBTO_placeholders)
        
        parameters = ab_id.copy()
        parameters.extend(idCL)
        parameters.extend(idBTO)
        exp_df = pd.read_sql(sql=idCL_BTO_query, params=parameters, con=conn)
    
    # print(exp_df)
    return exp_df