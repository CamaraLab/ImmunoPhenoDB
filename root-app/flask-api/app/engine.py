import warnings
import logging

import configparser
import mysql.connector # mysql-connector-python

import re

import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from statsmodels.regression.mixed_linear_model import MixedLMResults
from statsmodels.iolib.summary2 import Summary

import pickle
import networkx as nx
import requests

from pyspark.sql import SparkSession
from pyspark import SparkContext
import pymysql

# For STvEA
import networkx as nx
import numpy as np
from random import sample
import math

SCI_CRUNCH_BASE = "http://www.scicrunch.org"
SCI_RRID_ENDPOINT = "/resolver/"
SCI_FILE = ".json"
UNIPROT_BASE = "https://rest.uniprot.org"
UNIPROT_ENDPOINT = "/uniprotkb/search"
EBI_BASE = "https://www.ebi.ac.uk/ols/api/search"

def _config(filename: str = '/app/config/config.ini', 
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

#------------------------ Functions for get_tissues ------------------------#
def get_tissues() -> pd.DataFrame:
    """
    Retrieves a table of tissue types from a database

    Parameters:
        N/A
    
    Returns:
        tissues_df (pd.DataFrame): DataFrame containing tissue name IDs
            and tissue names
        
    """
    params = _config()
    conn = mysql.connector.connect(**params)

    tissues_query = """SELECT DISTINCT * FROM tissues;"""
    tissues_df = pd.read_sql(sql=tissues_query, con=conn)

    if conn is not None:
        conn.close()
        
    return tissues_df

#------------------------ Functions for get_antibodies ------------------------#
def get_antibodies_idCL(id_CLs: list,
                        idBTO: list = None,
                        idExperiment: list = None) -> list:
    """
    Retrieves a list of antibodies that have been used in cells
    categorized by certain cell types (idCLs), tissues (idBTO),
    and experiment (idExperiment)

    Parameters:
        id_CLs (list): cell types used to filter the
            search results of cells
        idBTO (list): tissue types used to filter the
            search results of cells
        idExperiment (list): database experiment IDs used to filter
            the search results of cells

    Returns:
        ab_list (list): antibodies used from the provided 
            parameters
    """
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()
    
    # Dynamically generate placeholders for parameterized query
    idCL_placeholders = ','.join(['%s'] * len(id_CLs))
    
    # No additional filtering
    if idBTO is None and idExperiment is None:
        ab_query = """SELECT DISTINCT antigen_expression.idAntibody
                    FROM antigen_expression
                    JOIN cells ON antigen_expression.idCell = cells.idCell
                    WHERE cells.idCL IN (%s);""" % (idCL_placeholders)
        
        cursor.execute(ab_query, (id_CLs))
    
    # Filtering based on tissue type (idBTO)
    elif idBTO is not None and idExperiment is None:
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
    
        ab_idBTO_query = """SELECT DISTINCT antigen_expression.idAntibody
                        FROM antigen_expression
                        JOIN cells ON antigen_expression.idCell = cells.idCell
                        JOIN experiments ON cells.idExperiment = experiments.idExperiment
                        WHERE cells.idCL IN (%s) AND experiments.idBTO IN (%s);""" % (idCL_placeholders, idBTO_placeholders)
        
        temp_idCLs = id_CLs.copy()
        temp_idCLs.extend(idBTO)

        cursor.execute(ab_idBTO_query, (temp_idCLs))
    
    # Filtering based on experiment id (idExperiment)
    elif idBTO is None and idExperiment is not None:
        idExperiment_placeholders = ','.join(['%s'] * len(idExperiment))
       
        ab_idExperiment_query = """SELECT DISTINCT antigen_expression.idAntibody
                                FROM antigen_expression
                                JOIN cells ON antigen_expression.idCell = cells.idCell
                                WHERE cells.idCL IN (%s) AND cells.idExperiment IN (%s);""" % (idCL_placeholders, idBTO_placeholders)

        temp_idCLs = id_CLs.copy()
        temp_idCLs.extend(idExperiment)
        cursor.execute(ab_idExperiment_query, (temp_idCLs))
    
    # Filtering based on both idBTO and idExperiment
    elif idBTO is not None and idExperiment is not None:
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
        idExperiment_placeholders = ','.join(['%s'] * len(idExperiment))
       
        ab_idBTO_idExp_query = """SELECT DISTINCT antigen_expression.idAntibody
                                FROM antigen_expression
                                JOIN cells ON antigen_expression.idCell = cells.idCell
                                JOIN experiments ON cells.idExperiment = experiments.idExperiment
                                WHERE cells.idCL IN (%s) AND cells.idExperiment IN (%s) AND experiments.idBTO IN (%s);""" % (idCL_placeholders, idExperiment_placeholders, idBTO_placeholders)
        
        temp_idCLs = id_CLs.copy()
        temp_idCLs.extend(idExperiment)
        temp_idCLs.extend(idBTO)

        cursor.execute(ab_idBTO_idExp_query, (temp_idCLs))
    
    res = cursor.fetchall()
    ab_list = [ab[0] for ab in res]
    
    if conn is not None:
        cursor.close()
        conn.close()
    
    return ab_list

def retrieve_lmm(summary_table: MixedLMResults) -> Summary:
    """
    Extracts the summarized mixed model regression results

    Parameters:
        summary_table (MixedLMResults): table generated after fitting
            a MixedLM model
    
    Returns:
        core (Summary): table containing coefficients, p values, etc 
    """
    # Return bottom half as a dataframe
    core = summary_table.summary().tables[1]

    return core

def spark_query_ab_lmm(antibody_id: str, node_fam_dict: dict) -> tuple:
    """
    Mapper function used when performing a Spark job to find LMM results,
    proportion of background cells for each cell type

    Finds appropriate cells given antibody and cell type restrictions
    and runs a linear mixed effects model. Also finds the proportion of
    cells expressed as background for each respective cell type

    Parameters:
        antibody_id (str): antibody ID
        node_fam_dict (dict): dictionary containing all descendants
            or family members for a given cell type (idCL). This dictionary
            is generated client side and sent to the server

            Example:
            Input (client): "CL:000236" 
            Output (client): 
                {'CL:0000236': ['CL:0000955',
                                'CL:0000788',
                                'CL:0000787',
                                'CL:0000816',
                                'CL:0001201',
                                'CL:0000785',
                                'CL:0000817']}
    
    Returns:
        antibody_id, core_copy, aggregate_df, bg_proportion (tuple):

        tuple containing:
            antibody_id (str): antibody ID
            core_copy (Summary): table containing results of LMM
            aggregate_df (pd.DataFrame): table containing background proportions
                divided by each cell type
            bg_proportion (float): proportion of cells classified as background
                belonging to cell types not specified by the user
    """
    # Extract all nodes in family of idCLs
    parents = list(node_fam_dict.keys())
    children = [item for sublist in list(node_fam_dict.values()) for item in sublist]
    family = parents + children

    # MySQL connection properties
    mysql_properties = _config()
    
    # Convert port from string to int in the dictionary
    mysql_properties['port'] = int(mysql_properties['port'])

    # Connect to MySQL
    connection = pymysql.connect(**mysql_properties)

    # Create a cursor object
    cursor = connection.cursor()

    # Get cells for LMM
    combined_query = """
        SELECT
            cells.idCell,
            cells.idCellOriginal,
            antigen_expression.normValue,
            cells.idExperiment,
            0 AS cellType -- marks cells that are associated with ab and cell type
        FROM
            antigen_expression
            INNER JOIN cells ON antigen_expression.idCell = cells.idCell
        WHERE
            antigen_expression.idAntibody = %s
            AND cells.idCL IN %s

        UNION ALL  

        SELECT
            cells.idCell,
            cells.idCellOriginal,
            antigen_expression.normValue,
            cells.idExperiment,
            1 AS cellType -- marks cells that are associated with ab but NOT cell type
        FROM
            antigen_expression
            INNER JOIN cells ON antigen_expression.idCell = cells.idCell
        WHERE
            antigen_expression.idAntibody = %s
            AND cells.idCL NOT IN %s;
    """

    # Convert cell_ids list to a tuple for both queries
    cell_ids_tuple = tuple(family)

    # Execute the combined query with parameters
    cursor.execute(combined_query, (antibody_id, cell_ids_tuple, antibody_id, cell_ids_tuple))

    # PART 1: get cells
    # Fetch the results into a Pandas DataFrame
    result_df = pd.DataFrame(cursor.fetchall(), columns=["idCell", "idCellOriginal", "normValue", "idExperiment", "cellType"])

    # EDGE CASE: If the provided idCLs cover all the cells that are used by that antibody
    # It will cause the LMM to fail
    # Skip this antibody by just returning a tuple() containing the antibody ID
    if len(set(result_df["cellType"])) <= 1:
        return (antibody_id, )
    
    # EDGE CASE: If all normalized values in result_df are background, it will cause the LMM to produce NaN values
    if (set(result_df["normValue"]) == {0}):
        return (antibody_id, )

    # PART 2: run LMM
    # Start LMM
    lmm = smf.mixedlm("normValue ~ cellType", result_df, groups=result_df["idExperiment"])
    summary_table = lmm.fit(method=["bfgs"])
    lmm_res = retrieve_lmm(summary_table)

    # Add in p-values
    p_values = summary_table.pvalues
    core_copy = lmm_res.copy(deep=True).drop(columns='P>|z|')
    core_copy['P>|z|'] = pd.Series(dtype="float64")
    core_copy.loc['Intercept', 'P>|z|'] = p_values[0]
    core_copy.loc['cellType', 'P>|z|'] = p_values[1]
    
    # PART 3: find background proportion of cells based on idCL
    with_query = """
        SELECT
            antigen_expression.idAntibody,
            cells.idCL,
            COUNT(antigen_expression.background) AS total_cells,
            SUM(antigen_expression.background = 1) AS bg_cells,
            SUM(antigen_expression.background = 1) / COUNT(antigen_expression.background) AS proportion
        FROM
            antigen_expression
            INNER JOIN cells ON antigen_expression.idCell = cells.idCell
            INNER JOIN (
                SELECT DISTINCT antigen_expression.idCell
                FROM antigen_expression
                WHERE antigen_expression.idAntibody = %s
            ) AS subquery_cells ON cells.idCell = subquery_cells.idCell
            INNER JOIN (
                SELECT DISTINCT cells.idCell
                FROM cells
                WHERE cells.idCell IN (
                    SELECT antigen_expression.idCell
                    FROM antigen_expression
                    WHERE antigen_expression.idAntibody = %s
                ) AND cells.idCL IN %s
            ) AS subquery_cl ON cells.idCell = subquery_cl.idCell
        WHERE
            antigen_expression.idAntibody = %s
        GROUP BY
            cells.idCL;
    """

    # Convert idCLs list to a tuple
    idCL_family_tuple = tuple(family)

    # Execute the query 
    cursor.execute(with_query, (antibody_id, antibody_id, idCL_family_tuple, antibody_id))
    
    prop_cells_with_idCL_df = pd.DataFrame(cursor.fetchall(), 
                                           columns=["idAntibody", 
                                                    "idCL", 
                                                    "total_cells", 
                                                    "bg_cells", 
                                                    "proportion"])
    
    new_df = prop_cells_with_idCL_df.set_index('idCL')
    df_idCLs = list(prop_cells_with_idCL_df['idCL'])
    aggregate_df = pd.DataFrame(data=0, index=node_fam_dict.keys(), columns=['total_cells', 'bg_cells'])

    # We want to aggregate idCLs that are parent-descendant related
    for key, value in node_fam_dict.items():
        # If our parent key has a value in the df
        if key in df_idCLs:
            # Add their total cells 
            aggregate_df.loc[key]['total_cells'] += int(new_df.loc[key]['total_cells'])
            # Add their background cells
            aggregate_df.loc[key]['bg_cells'] += int(new_df.loc[key]['bg_cells'])
            
        for child in value:
            if child in df_idCLs:
                # Add their total cells to their parent's total_cells column
                aggregate_df.loc[key]['total_cells'] += int(new_df.loc[child]['total_cells'])
                # Add their background cells
                aggregate_df.loc[key]['bg_cells'] += int(new_df.loc[child]['bg_cells'])
                
        # WARNING: if there are 0 cells for a given idCL for an antibody
        if aggregate_df.loc[key]['total_cells'] == 0:
            logging.warning(f"No cells found for type {key} in {antibody_id}. "
                           "Background percentage will be set to NaN. ")

    # Calculate the proportions of background cells into a new column
    aggregate_df['proportion'] = aggregate_df['bg_cells']/aggregate_df['total_cells']

    without_query = """
        SELECT
            COUNT(*) AS total,
            SUM(antigen_expression.background = 1) AS bg_cells,
            SUM(antigen_expression.background = 1) / COUNT(*) AS proportion
        FROM
            antigen_expression
            INNER JOIN cells ON antigen_expression.idCell = cells.idCell
        WHERE
            antigen_expression.idCell IN (
                SELECT antigen_expression.idCell
                FROM antigen_expression
                INNER JOIN cells ON antigen_expression.idCell = cells.idCell
                WHERE
                    antigen_expression.idAntibody = %s
                    AND cells.idCL NOT IN %s
            )
            AND antigen_expression.idAntibody = %s;
    """

    # Execute the query 
    cursor.execute(without_query, (antibody_id, idCL_family_tuple, antibody_id))
    prop_cells_without_idCL_df = pd.DataFrame(cursor.fetchall(), columns=["total", "bg_cells", "proportion"])
    bg_proportion = float(prop_cells_without_idCL_df['proportion'].values[0])
    
    # Close the cursor and MySQL connection
    cursor.close()
    connection.close()

    return (antibody_id, core_copy, aggregate_df, bg_proportion)

def run_spark_query_ab_lmm(antibody_ids: list, node_fam_dict: dict) -> dict:
    """
    Performs a Spark job for a list of antibodies. Collects all results and
    reduces the results into a dictionary. 

    Spark job will run spark_query_ab_lmm in parallel for antibody_ids. 

    Parameters:
        antibody_ids (list): list of antibody IDs
        node_fam_dict (dict): dictionary containing all descendants
            or family members for a given cell type (idCL). This dictionary
            is generated client side and sent to the server

            Example:
            Input (client): "CL:000236" 
            Output (client): 
                {'CL:0000236': ['CL:0000955',
                                'CL:0000788',
                                'CL:0000787',
                                'CL:0000816',
                                'CL:0001201',
                                'CL:0000785',
                                'CL:0000817']}

    Returns:
        ab_lmm_dict (dict): dictionary containing results of Spark job
        Example: 
            {'AB_2800725': (LMM, aggregate_df, bg_proportion)}
            key: antibody ID
            value: tuple containing:
                LMM (Summary): table containing results of LMM
                aggregate_df (pd.DataFrame): table containing background proportions
                    divided by each cell type
                bg_proportion (float): proportion of cells classified as background
                    belonging to cell types not specified by the user
    """
    # Create a Spark session
    spark = SparkSession.builder.appName("ParallelMySQLQuery").getOrCreate()

    # We've already provided the idCLs to iterate over
    num_parallel_queries = 100

    # Create an RDD 
    rdd_abs = spark.sparkContext.parallelize(antibody_ids, num_parallel_queries)

    # Run the spark_query_ab_lmm in parallel and store in a list of dataframes
    results = rdd_abs.map(lambda x: spark_query_ab_lmm(x, node_fam_dict)).collect()

    # Stop spark
    spark.stop()

    # Convert this list of results into a dictionary
    # where keys: antibody id, value: (LMM, aggregate_df, not_idCL_proportion)
    # Organize results based on antibody
    ab_lmm_dict = {}
    for result in results:
        # EDGE CASE: If there were no results for a given antibody, 
        # it will be a tuple with only the antibody ID
        # Skip it
        antibody_name = result[0]
        if len(result) == 1:
            continue
            
        if antibody_name in antibody_ids:
            ab_lmm_dict[antibody_name] = (result[1], result[2], result[3])

    return ab_lmm_dict

def lmm_ab_list_to_df(lmm_results: dict, parents: list) -> pd.DataFrame:
    """
    Converts output from Spark job into the final dataframe format
    
    Parameters:
        lmm_results (dict): results from Spark job
            Example: 
            {'AB_2800725': (LMM, aggregate_df, bg_proportion)}
            key: antibody ID
            value: tuple containing:
                LMM (Summary): table containing results of LMM
                aggregate_df (pd.DataFrame): table containing background proportions
                    divided by each cell type
                bg_proportion (float): proportion of cells classified as background
                    belonging to cell types not specified by the user
            parents (list): list of cell types (idCL) without any descendants

    Returns:
        final_df (pd.DataFrame): dataframe containing the final output of the Spark job
            containing LMM results and background proportions for a given antibody and
            cell type
    """
    cellType_coeff = []
    stdError = []
    p_values = []

    antibodies = list(lmm_results.keys())
    
    idCL_dict = {k: [] for k in parents}
    idCL_dict['other'] = []
    
    # Loop over each antibody_id (key) in the outer dictionary
    for antibody_id, results in lmm_results.items():
        # results is: LMM table, aggregate_df, proportion of background cells not in idCLs list
        lmm_table = results[0]
        aggregate_df = results[1]
        bg_proportion = results[2]
        
        # Taking the negative of the coefficient
        cellType_coeff.append(-1 * float(lmm_table.loc['cellType']['Coef.']))
        stdError.append(float(lmm_table.loc['cellType']['Std.Err.']))
        p_values.append(float(lmm_table.loc['cellType']['P>|z|']))

        for idCL in parents:
            # Take (1 - proportion_background), which will instead show proportion of signal cells
            proportion_background = aggregate_df.loc[idCL]['proportion']
            idCL_dict[idCL].append(1 - proportion_background)

        idCL_dict['other'].append(bg_proportion)
            
    # Apply Benjamini Hochberg correction to all the p_values
    q_values = multipletests(p_values, method="fdr_bh")[1]

    table_dict = {'coeff': cellType_coeff,
            'stderr': stdError, 
            'p_val': p_values,
            'q_val': q_values}

    entire_dict = {**table_dict, **idCL_dict}
    
    final_df = pd.DataFrame(index=antibodies, data=entire_dict)

    # Add a readable column at the end by doing a query to the database
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()
    
    ab_query = """SELECT antibodies.abTarget
                FROM antibodies
                WHERE antibodies.idAntibody = (%s);"""
    
    # Add cellType column of names from id
    ab_conv = []
    for ab in antibodies:
        cursor.execute(ab_query, (ab,))
        ab_name = cursor.fetchone()[0]
        ab_conv.append(ab_name)
    
    final_df.insert(0, 'target', ab_conv)
    final_df.sort_values('coeff', ascending=False, inplace=True)

    return final_df

def spark_query_custom_ab_lmm(antibody_id: str, 
                              node_fam_dict: dict, 
                              custom_background_fam_dict: dict):
    # Extract all nodes in family, used as our target input
    target_parents = list(node_fam_dict.keys())
    target_children = [item for sublist in list(node_fam_dict.values()) for item in sublist]
    target_family = target_parents + target_children

    # Extract all nodes in the custom background cell types
    background_parents = list(custom_background_fam_dict.keys())
    background_children = [item for sublist in list(custom_background_fam_dict.values()) for item in sublist]
    background_family = background_parents + background_children

    # MySQL connection properties
    mysql_properties = _config()
    
    # Convert port from string to int in the dictionary
    mysql_properties['port'] = int(mysql_properties['port'])

    # Connect to MySQL
    connection = pymysql.connect(**mysql_properties)

    # Create a cursor object
    cursor = connection.cursor()

    # Get cells for LMM
    combined_query = """
        SELECT
            cells.idCell,
            cells.idCellOriginal,
            antigen_expression.normValue,
            cells.idExperiment,
            0 AS cellType -- marks cells that are associated with ab and cell type
        FROM
            antigen_expression
            INNER JOIN cells ON antigen_expression.idCell = cells.idCell
        WHERE
            antigen_expression.idAntibody = %s
            AND cells.idCL IN %s

        UNION ALL  -- Use UNION ALL to combine the results of both queries

        SELECT
            cells.idCell,
            cells.idCellOriginal,
            antigen_expression.normValue,
            cells.idExperiment,
            1 AS cellType -- marks background cells based on provided cell type
        FROM
            antigen_expression
            INNER JOIN cells ON antigen_expression.idCell = cells.idCell
        WHERE
            antigen_expression.idAntibody = %s
            AND cells.idCL IN %s;
    """

    # Convert target idCLs to tuple
    target_family_tuple = tuple(target_family)

    # Convert background idCLs to tuple
    background_family_tuple = tuple(background_family)

    # Execute the combined query with parameters
    cursor.execute(combined_query, (antibody_id, target_family_tuple, antibody_id, background_family_tuple))

    # Fetch the results into a Pandas DataFrame
    result_df = pd.DataFrame(cursor.fetchall(), columns=["idCell", "idCellOriginal", 
                                                         "normValue", "idExperiment", "cellType"])
    
    # EDGE CASE: If the provided idCLs cover all the cells that are used by that antibody, it will cause the LMM to fail
    # Skip this antibody by just returning an empty tuple()
    if len(set(result_df["cellType"])) <= 1:
        return (antibody_id, )
    
    # EDGE CASE: If all normalized values in result_df are background, it will cause the LMM to produce NaN values
    if (set(result_df["normValue"]) == {0}):
        return (antibody_id, )

    # # Start LMM
    lmm = smf.mixedlm("normValue ~ cellType", result_df, groups=result_df["idExperiment"])
    summary_table = lmm.fit(method=["bfgs"])
    lmm_res = retrieve_lmm(summary_table)

    # Add in p-values manually (with increased precision)
    p_values = summary_table.pvalues
    core_copy = lmm_res.copy(deep=True).drop(columns='P>|z|')
    core_copy['P>|z|'] = pd.Series(dtype="float64")
    core_copy.loc['Intercept', 'P>|z|'] = p_values[0]
    core_copy.loc['cellType', 'P>|z|'] = p_values[1]
    
    # We now want to find the percentage of background cells for our target cell types
    with_query = """
        SELECT
            antigen_expression.idAntibody,
            cells.idCL,
            COUNT(antigen_expression.background) AS total_cells,
            SUM(antigen_expression.background = 1) AS bg_cells,
            SUM(antigen_expression.background = 1) / COUNT(antigen_expression.background) AS proportion
        FROM
            antigen_expression
            INNER JOIN cells ON antigen_expression.idCell = cells.idCell
            INNER JOIN (
                SELECT DISTINCT antigen_expression.idCell
                FROM antigen_expression
                WHERE antigen_expression.idAntibody = %s
            ) AS subquery_cells ON cells.idCell = subquery_cells.idCell
            INNER JOIN (
                SELECT DISTINCT cells.idCell
                FROM cells
                WHERE cells.idCell IN (
                    SELECT antigen_expression.idCell
                    FROM antigen_expression
                    WHERE antigen_expression.idAntibody = %s
                ) AND cells.idCL IN %s
            ) AS subquery_cl ON cells.idCell = subquery_cl.idCell
        WHERE
            antigen_expression.idAntibody = %s
        GROUP BY
            cells.idCL;
    """

    # Execute the query 
    cursor.execute(with_query, (antibody_id, antibody_id, target_family_tuple, antibody_id))
    
    prop_cells_with_idCL_df = pd.DataFrame(cursor.fetchall(), 
                                           columns=["idAntibody", "idCL", "total_cells", "bg_cells", "proportion"])
    
    new_df = prop_cells_with_idCL_df.set_index('idCL')
    df_idCLs = list(prop_cells_with_idCL_df['idCL'])
    target_aggregate_df = pd.DataFrame(data=0, index=node_fam_dict.keys(), columns=['total_cells', 'bg_cells'])

    # We want to aggregate idCLs that are parent-descendant related
    for key, value in node_fam_dict.items():
        # If our parent key has a value in the df
        if key in df_idCLs:
            # Add their total cells 
            target_aggregate_df.loc[key]['total_cells'] += int(new_df.loc[key]['total_cells'])
            # Add their background cells
            target_aggregate_df.loc[key]['bg_cells'] += int(new_df.loc[key]['bg_cells'])
            
        for child in value:
            if child in df_idCLs:
                # Add their total cells to their parent's total_cells column
                target_aggregate_df.loc[key]['total_cells'] += int(new_df.loc[child]['total_cells'])
                # Add their background cells
                target_aggregate_df.loc[key]['bg_cells'] += int(new_df.loc[child]['bg_cells'])
                
        # WARNING: if there are 0 cells for a given idCL for an antibody
        if target_aggregate_df.loc[key]['total_cells'] == 0:
            logging.warning(f"No cells found for type {key} in {antibody_id}. "
                           "Background percentage will be set to NaN. ")

    # Calculate the proportions of background cells for our target idCL into a new column
    target_aggregate_df['proportion'] = target_aggregate_df['bg_cells']/target_aggregate_df['total_cells']

    ########################################################################################
    # We have to find the proportion of background cells for each of the custom background idCLs
    without_query = """
        SELECT
            antigen_expression.idAntibody,
            cells.idCL,
            COUNT(antigen_expression.background) AS total_cells,
            SUM(antigen_expression.background = 1) AS bg_cells,
            SUM(antigen_expression.background = 1) / COUNT(antigen_expression.background) AS proportion
        FROM
            antigen_expression
            INNER JOIN cells ON antigen_expression.idCell = cells.idCell
            INNER JOIN (
                SELECT DISTINCT antigen_expression.idCell
                FROM antigen_expression
                WHERE antigen_expression.idAntibody = %s
            ) AS subquery_cells ON cells.idCell = subquery_cells.idCell
            INNER JOIN (
                SELECT DISTINCT cells.idCell
                FROM cells
                WHERE cells.idCell IN (
                    SELECT antigen_expression.idCell
                    FROM antigen_expression
                    WHERE antigen_expression.idAntibody = %s
                ) AND cells.idCL IN %s
            ) AS subquery_cl ON cells.idCell = subquery_cl.idCell
        WHERE
            antigen_expression.idAntibody = %s
        GROUP BY
            cells.idCL;
    """

    # Execute the query 
    cursor.execute(without_query, (antibody_id, antibody_id, background_family_tuple, antibody_id))
    
    prop_cells_without_idCL_df = pd.DataFrame(cursor.fetchall(), 
                                           columns=["idAntibody", "idCL", "total_cells", "bg_cells", "proportion"])
    
    without_new_df = prop_cells_without_idCL_df.set_index('idCL')

    # We also need to aggregate this now too
    background_df_idCLs = list(prop_cells_without_idCL_df['idCL'])
    background_aggregate_df = pd.DataFrame(data=0, 
                                           index=custom_background_fam_dict.keys(), 
                                           columns=['total_cells', 'bg_cells'])

    # We want to aggregate idCLs that are parent-descendant related
    for key, value in custom_background_fam_dict.items():
        # If our parent key has a value in the df
        if key in background_df_idCLs:
            # Add their total cells 
            background_aggregate_df.loc[key]['total_cells'] += int(without_new_df.loc[key]['total_cells'])
            # Add their background cells
            background_aggregate_df.loc[key]['bg_cells'] += int(without_new_df.loc[key]['bg_cells'])
            
        for child in value:
            if child in background_df_idCLs:
                # Add their total cells to their parent's total_cells column
                background_aggregate_df.loc[key]['total_cells'] += int(without_new_df.loc[child]['total_cells'])
                # Add their background cells
                background_aggregate_df.loc[key]['bg_cells'] += int(without_new_df.loc[child]['bg_cells'])
                
        # WARNING: if there are 0 cells for a given idCL for an antibody
        if background_aggregate_df.loc[key]['total_cells'] == 0:
            logging.warning(f"No cells found for type {key} in {antibody_id}. "
                           "Background percentage will be set to NaN. ")

    # Calculate the proportions of background cells into a new column
    background_aggregate_df['proportion'] = background_aggregate_df['bg_cells']/background_aggregate_df['total_cells']
    
    # Close the cursor and MySQL connection
    cursor.close()
    connection.close()

    return (antibody_id, core_copy, target_aggregate_df, background_aggregate_df)

def run_spark_query_custom_ab_lmm(antibody_ids: list, 
                                  node_fam_dict: dict, 
                                  custom_background_fam_dict: dict):
    # Create a Spark session
    spark = SparkSession.builder.appName("ParallelMySQLQuery").getOrCreate()

    # We've already provided the idCLs to iterate over
    num_parallel_queries = 100

    # Create an RDD 
    rdd_abs = spark.sparkContext.parallelize(antibody_ids, num_parallel_queries)

    # Run the spark_query_ab_lmm in parallel and store in a list of dataframes
    results = rdd_abs.map(lambda x: spark_query_custom_ab_lmm(x, node_fam_dict, 
                                                              custom_background_fam_dict)).collect()
    
    # Stop spark
    spark.stop()

    # Convert this list of results into a dictionary
    # where keys: antibody id, value: (LMM, target_aggregate_df, background_aggregate_df)
    # Organize results based on antibody
    ab_lmm_dict = {}
    for result in results:
        # EDGE CASE: If there were no results for a given antibody, it will be a tuple with only the antibody ID
        # Skip it
        antibody_name = result[0]
        if len(result) == 1:
            continue
            
        if antibody_name in antibody_ids:
            ab_lmm_dict[antibody_name] = (result[1], result[2], result[3])

    return ab_lmm_dict

def lmm_custom_ab_list_to_df(lmm_results: dict, 
                             target_idCL_parents: list, 
                             background_idCL_parents: list) -> pd.DataFrame:
    
    cellType_coeff = []
    stdError = []
    p_values = []

    antibodies = list(lmm_results.keys())
    
    target_idCL_dict = {k: [] for k in target_idCL_parents}
    background_idCL_dict = {k: [] for k in background_idCL_parents}
    
    # Loop over each antibody_id (key) in the outer dictionary
    for antibody_id, results in lmm_results.items():
        # results is: LMM table, target_aggregate_df, background_aggregate_df
        lmm_table = results[0]
        target_aggregate_df = results[1]
        background_aggregate_df = results[2]

        # Add LMM results in
        # Taking the negative of the coefficient
        cellType_coeff.append(-1 * float(lmm_table.loc['cellType']['Coef.']))
        stdError.append(float(lmm_table.loc['cellType']['Std.Err.']))
        p_values.append(float(lmm_table.loc['cellType']['P>|z|']))

        # Add target proportions in
        for idCL in target_idCL_parents:
            # Take (1 - proportion_background), which will instead show proportion of signal cells
            proportion = target_aggregate_df.loc[idCL]['proportion']
            target_idCL_dict[idCL].append(1 - proportion)

        # Add background proportions in 
        for idCL in background_idCL_parents:
            # Take (1 - proportion_background), which will instead show proportion of signal cells
            proportion = background_aggregate_df.loc[idCL]['proportion']
            background_idCL_dict[idCL].append(1 - proportion)
            
    # Apply Benjamini Hochberg correction to all the p_values
    q_values = multipletests(p_values, method="fdr_bh")[1]

    # Combining all parts together into a pandas dataframe
    table_dict = {'coeff': cellType_coeff,
            'stderr': stdError, 
            'p_val': p_values,
            'q_val': q_values}
    
    entire_dict = {**table_dict, **target_idCL_dict, **background_idCL_dict}
    final_df = pd.DataFrame(index=antibodies, data=entire_dict)

    # Add a readable column at the end by doing a query to the database
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()
    
    ab_query = """SELECT antibodies.abTarget
                FROM antibodies
                WHERE antibodies.idAntibody = (%s);"""
    
    # Add cellType column of names from id
    ab_conv = []
    for ab in antibodies:
        cursor.execute(ab_query, (ab,))
        ab_name = cursor.fetchone()[0]
        ab_conv.append(ab_name)
    
    final_df.insert(0, 'target', ab_conv)

    # Change sort from ascending to descending (since we took -1 * coefficient)
    final_df.sort_values('coeff', ascending=False, inplace=True)

    return final_df

def get_antibodies(node_fam_dict: dict,
                   custom_background_fam_dict: dict = None, 
                   idBTO: list = None,
                   idExperiment: list = None) -> pd.DataFrame:
    """
    Wrapper function to call the required queries and Spark job to find antibodies
    that mark a provided list of cell types

    Parameters:
        node_fam_dict (dict): dictionary containing all descendants
            or family members for a given cell type (idCL). This dictionary
            is generated client side and sent to the server

            Example:
            Input (client): "CL:000236" 
            Output (client): 
                {'CL:0000236': ['CL:0000955',
                                'CL:0000788',
                                'CL:0000787',
                                'CL:0000816',
                                'CL:0001201',
                                'CL:0000785',
                                'CL:0000817']}
        idBTO (list): list of tissue type IDs used for restricting the search query
        idExperiment (list): list of database experiment IDs used for restricting
            search query
    
    Returns:
        result_df (pd.DataFrame): formatted dataframe containing LMM results for
            each antibody
    """
    # Extract all nodes in family, used as our target input
    target_parents = list(node_fam_dict.keys())
    target_children = [item for sublist in list(node_fam_dict.values()) for item in sublist]
    target_family = target_parents + target_children

    # Retrieve unique antibodies from idCLs
    antibodies = get_antibodies_idCL(target_family, idBTO, idExperiment)

    # If background idCLs are provided, call custom query
    if custom_background_fam_dict is not None:
        # Extract all nodes in the custom background cell types
        background_parents = list(custom_background_fam_dict.keys())
        background_children = [item for sublist in list(custom_background_fam_dict.values()) for item in sublist]
        background_family = background_parents + background_children

         # Run custom Spark job on antibodies
        ab_lmm_output = run_spark_query_custom_ab_lmm(antibodies, node_fam_dict, custom_background_fam_dict)

    else:
        # Run regular Spark job on antibodies
        ab_lmm_output = run_spark_query_ab_lmm(antibodies, node_fam_dict)

    if len(ab_lmm_output) == 0:
        return ("No results found. Please retry with fewer cell types or different tissue and experiment parameters.")

    # Likewise for final output
    if custom_background_fam_dict is not None:
        # Convert output into a readable dataframe
        result_df = lmm_custom_ab_list_to_df(ab_lmm_output, target_parents, background_parents)
    else:
        result_df = lmm_ab_list_to_df(ab_lmm_output, target_parents) 

    return result_df

#------------------------ Functions for plot_antibodies ------------------------#
def group_normvalue(field:str, plot_values: pd.DataFrame) -> pd.DataFrame:
    """
    Helper function to aggregate normalized values based on column value

    Parameters:
        field (str): value in column to group by
        plot_values (pd.DataFrame): dataframe containing normalized values for 
            each cell with an antibody or cell type (idCL)
    
    Returns:
        grouped_df (pd.DataFrame): dataframe containing a column with
            antibody IDs, and a column holding lists of normalized values 
            from that antibody id or cell type (idCL)
    """
    grouped_df = plot_values.groupby(f'{field}')['normValue'].agg(list).reset_index()
    return grouped_df

def summary_statistics(field:str, grouped_df:pd.DataFrame) -> pd.DataFrame:
    """
    Calculates descriptive statistics for a list of normalized values

    Parameters:
        field (str): can be either "idAntibody" or "idCL"
        grouped_df (pd.DataFrame): dataframe containing normalized values,
            grouped by either antibody ID or cell type
    
    Returns:
        df_summary (pd.DataFrame): dataframe containing: 
            mean, quartile 1, median, quartile 3, min, and max values
    """
    # Calculate summary statistics for each "idAntibody"
    summary_stats = []
    for index, row in grouped_df.iterrows():
        id_field = row[f'{field}']
        norm_values = row['normValue']
        
        # Calculate summary statistics (mean, quartiles, etc.)
        mean_value = pd.Series(norm_values).mean()
        q1 = pd.Series(norm_values).quantile(0.25)
        median = pd.Series(norm_values).median()
        q3 = pd.Series(norm_values).quantile(0.75)
        min_value = pd.Series(norm_values).min()
        max_value = pd.Series(norm_values).max()
        
        summary_stats.append({f'{field}': id_field, 
                              'mean': mean_value, 
                              'q1': q1, 
                              'median': median, 
                              'q3': q3,
                              'min': min_value, 
                              'max': max_value})
    
    # Create a DataFrame from the summary statistics
    df_summary = pd.DataFrame(summary_stats)

    return df_summary

def plot_antibodies(ab_ids: list, 
                    id_CLs: list)-> pd.DataFrame:
    """
    Queries a database to find normalized values based on cells found in
    a certain subset of antibodies and cell types

    Processes queried data into descriptive statistics used for plotting
    on the client.

    Parameters:
        ab_ids (list): list of antibody IDs used for querying
        id_CLs (list): list of cell type IDs used for querying
    
    Returns:
        summary_stats (pd.DataFrame): dataframe containing the min, max, 
            mean, etc for a particular antibody's normalized values
    """
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)
  
    # Dynamically generate placeholders for parameterized query
    ab_placeholders = ','.join(['%s'] * len(ab_ids))
    idCL_placeholders = ','.join(['%s'] * len(id_CLs)) # find descendants on client first, then send them here
    
    with_query = """SELECT cells.idCell, 
                        cells.idCellOriginal,
                        antigen_expression.normValue,
                        cells.idExperiment,
                        antigen_expression.idAntibody,
                        cells.idCL
                    FROM cells
                    INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                    WHERE cells.idCell IN (
                        SELECT antigen_expression.idCell
                        FROM antigen_expression 
                        WHERE antigen_expression.idAntibody IN (%s)
                    ) AND cells.idCL IN (%s)
                        AND antigen_expression.idAntibody IN (%s);"""  % (ab_placeholders, idCL_placeholders, ab_placeholders)

    parameters = ab_ids.copy()
    parameters.extend(id_CLs)
    parameters.extend(ab_ids)
    
    cells_with_idCL_df = pd.read_sql(sql=with_query, params=parameters, con=conn)

    grouped_normvalue_by_antibody = group_normvalue("idAntibody", cells_with_idCL_df)

    summary_stats = summary_statistics("idAntibody", grouped_normvalue_by_antibody)
    
    if conn is not None:
        conn.close()
    
    return summary_stats

#------------------------ Functions for get_celltype ------------------------#
def get_idCLs_ab_ids(ab_id: str,
                     idBTO: list = None, 
                     idExperiment: list = None) -> list:
    """
    Retrieves a list of cell types that have been used in cells
    categorized by a certain antibody (ab_id), tissues (idBTO),
    and experiment (idExperiment)

    Parameters:
        ab_id (str): a single antibody ID of interest to narrow
            cells in search result
        idBTO (list): tissue types used to filter the
            search results of cells
        idExperiment (list): database experiment IDs used to filter
            the search results of cells
    
    Returns:
        idCL_list (list): list of cell type IDs
    """
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()
    
    # No additional filtering
    if idBTO is None and idExperiment is None:
        idCL_query = """SELECT DISTINCT cells.idCL
                    FROM cells
                    INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                    WHERE antigen_expression.idAntibody = (%s);""" % ('%s')

        cursor.execute(idCL_query, (ab_id, ))
    
    # Filtering based on idBTO
    elif idBTO is not None and idExperiment is None:
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))

        idCL_idBTO_query = """SELECT DISTINCT cells.idCL
                            FROM cells
                            INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
                            INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                            WHERE experiments.idBTO IN (%s) AND antigen_expression.idAntibody = (%s);""" % (idBTO_placeholders, '%s')

        temp_ab = idBTO.copy()
        temp_ab.append(ab_id)
        
        cursor.execute(idCL_idBTO_query, (temp_ab))
    
    # Filtering based on idExperiment
    elif idBTO is None and idExperiment is not None:
        idExperiment_placeholders = ','.join(['%s'] * len(idExperiment))
        
        idCL_idExp_query = """SELECT DISTINCT cells.idCL
                            FROM cells
                            INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                            WHERE cells.idExperiment IN (%s) AND antigen_expression.idAntibody = (%s);""" % (idExperiment_placeholders, '%s')

        temp_ab = idExperiment.copy()
        temp_ab.append(ab_id)
        
        cursor.execute(idCL_idExp_query, (temp_ab))
    
    # Filtering based on idBTO and idExperiment
    elif idBTO is not None and idExperiment is not None:
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
        idExperiment_placeholders = ','.join(['%s'] * len(idExperiment))

        idCL_idBTO_idExp_query = """SELECT DISTINCT cells.idCL
                                FROM cells
                                INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
                                INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                                WHERE experiments.idBTO IN (%s) AND cells.idExperiment IN (%s) AND antigen_expression.idAntibody = (%s);"""
        
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

def spark_query_cell_type_lmm(antibody_id: str, idCL: str) -> tuple:
    """
    Mapper function used when performing the second Spark job 
    in the pipeline for get_celltypes

    Finds appropriate cells given an antibody and cell type 
    and runs a linear mixed effects model. Also finds the proportion of
    cells expressed as background for the specified cell type

    Parameters:
        antibody_id (str): antibody ID used to filter results
        idCL (str): cell type ID used to filter results
    
    Returns:
        antibody_id, idCL, core_copy, proportion_background (tuple):

        tuple containing:
            antibody_id (str): antibody ID
            idCL (str): cell type ID
            core_copy (Summary): table containing results of LMM
            bg_proportion (float): proportion of cells classified as background
                belonging to the specific cell type and antibody
    """
    # MySQL connection properties
    mysql_properties = _config()
    
    # Convert port from string to int in the dictionary
    mysql_properties['port'] = int(mysql_properties['port'])

    # Connect to MySQL
    connection = pymysql.connect(**mysql_properties)

    # Create a cursor object
    cursor = connection.cursor()

    # Combined MySQL query with placeholders
    combined_query = """
        SELECT
            cells.idCell,
            cells.idCellOriginal,
            antigen_expression.normValue,
            cells.idExperiment,
            antigen_expression.idAntibody,
            0 AS cellType -- marks cells that are associated with ab and cell type
        FROM
            antigen_expression
        INNER JOIN cells ON antigen_expression.idCell = cells.idCell
        WHERE
            antigen_expression.idAntibody = %s
            AND cells.idCell IN (
                SELECT
                    antigen_expression.idCell
                FROM
                    antigen_expression
                INNER JOIN cells ON antigen_expression.idCell = cells.idCell
                WHERE
                    antigen_expression.idAntibody = %s
                    AND cells.idCL = %s
            )

        UNION ALL
        
        SELECT
            cells.idCell,
            cells.idCellOriginal,
            antigen_expression.normValue,
            cells.idExperiment,
            antigen_expression.idAntibody,
            1 AS cellType -- marks cells that are associated with ab but NOT cell type
        FROM
            antigen_expression
        INNER JOIN cells ON antigen_expression.idCell = cells.idCell
        WHERE
            antigen_expression.idAntibody = %s
            AND cells.idCell IN (
                SELECT
                    antigen_expression.idCell
                FROM
                    antigen_expression
                INNER JOIN cells ON antigen_expression.idCell = cells.idCell
                WHERE
                    antigen_expression.idAntibody = %s
                    AND cells.idCL != %s
            );
            """

    # Execute the combined query with parameters
    cursor.execute(combined_query, (antibody_id, antibody_id, idCL, antibody_id, antibody_id, idCL))

    # PART 1: get cells
    # Fetch the results into a Pandas DataFrame
    result_df = pd.DataFrame(cursor.fetchall(), columns=["idCell", "idCellOriginal", "normValue", "idExperiment", "idAntibody", "cellType"])
    
    # PART 2: run LMM
    # Start LMM
    lmm = smf.mixedlm("normValue ~ cellType", result_df, groups=result_df["idExperiment"])
    summary_table = lmm.fit(method=["bfgs"])
    lmm_res = retrieve_lmm(summary_table)

    # Add in p-values 
    p_values = summary_table.pvalues
    core_copy = lmm_res.copy(deep=True).drop(columns='P>|z|')
    core_copy['P>|z|'] = pd.Series(dtype="float64")
    core_copy.loc['Intercept', 'P>|z|'] = p_values[0]
    core_copy.loc['cellType', 'P>|z|'] = p_values[1]

    # PART 3: find proportion of background cells for this exact antibody and cell type
    background_query = """
        SELECT
            antigen_expression.idAntibody,
            cells.idCL,
            COUNT(antigen_expression.background) as total_cells,
            SUM(IF(antigen_expression.background = 1, 1, 0)) as bg_cells,
            IFNULL(SUM(IF(antigen_expression.background = 1, 1, 0)) / NULLIF(COUNT(antigen_expression.background), 0), 0) as proportion
        FROM
            antigen_expression
        INNER JOIN
            cells ON antigen_expression.idCell = cells.idCell
        WHERE
            antigen_expression.idAntibody = (%s)
            AND cells.idCell IN (
                SELECT
                    cells.idCell
                FROM
                    cells
                WHERE
                    cells.idCell IN (
                        SELECT
                            antigen_expression.idCell
                        FROM
                            antigen_expression
                        WHERE
                            antigen_expression.idAntibody = (%s)
                    )
                    AND cells.idCL = (%s)
            )
        GROUP BY
            antigen_expression.idAntibody, cells.idCL;"""
    
    cursor.execute(background_query, (antibody_id, antibody_id, idCL))
    # Fetch the results into a Pandas DataFrame
    bg_df = pd.DataFrame(cursor.fetchall(), columns=["idAntibody", "idCL", "total_cells", "bg_cells", "proportion"])
    bg_df2 = bg_df.set_index("idAntibody")

    # Better accuracy is to compute it ourselves (MySQL does rounding)
    if bg_df2.loc[antibody_id]["total_cells"] == 0:
        logging.warning(f"No cells found for type {idCL} in {antibody_id}. "
                       "Proportion will be set to NaN")

    proportion_background = float(bg_df2["bg_cells"]/bg_df2["total_cells"])    
    
    # Close the cursor and MySQL connection
    cursor.close()
    connection.close()

    return (antibody_id, idCL, core_copy, proportion_background)

def run_spark_cell_type_from_ab(antibody_ids:list,
                                idBTO:list = None,
                                idExperiment:list = None) -> list:
    """
    Performs the first Spark job in the get_celltypes pipeline.
    For a list of antibodies, collects all results from the Spark job and
    reduces the results into a list.

    Spark job will run get_idCLs_ab_ids in parallel, which finds all cell types
    of cells that have been used with the antibodies provided in the parameter,
    optionally restricted by tissue type and experiment ID.

    Parameters:
        antibody_ids (list): antibody IDs used to find cell types
        idBTO (list): tissue type IDs to restrict cell type search results
        idExperiment: experiment IDs to restrict cell type search results

    Returns:
        pairs_list (list): list containing tuples of each antibody with a cell type
            Example:
                [(antibody_id, idCL), (antibody_id, idCL) ... ]
    """
    # Create a Spark session
    spark = SparkSession.builder.appName("ParallelMySQLQuery").getOrCreate()

    # We're already provided the idCLs to iterate over
    num_parallel_queries = 10

    # Create an RDD 
    rdd_abs = spark.sparkContext.parallelize(antibody_ids, num_parallel_queries)

    # Define a new map function to include our additional parameters
    def map_function(ab_id):
        return get_idCLs_ab_ids(ab_id, idBTO=idBTO, idExperiment=idExperiment)
    
    results = rdd_abs.map(map_function).collect()

    # Stop spark
    spark.stop()

    old_list = []
    pairs_list = []

    for idx, result_list in enumerate(results):
        ab_name = antibody_ids[idx]
        old_list.extend(result_list)
        for idCL in result_list:
            pairs_list.append((ab_name, idCL))

    return pairs_list

def run_spark_cell_type_lmm(antibody_ids: list, pairs_list: list) -> dict:
    """
    Performs the second Spark job in the get_celltypes pipeline using 
    the output from the first Spark job for the pairs_list parameter. 

    Spark job will run run_spark_cell_type_from_ab in parallel, which
    finds appropriate cells given an antibody and cell type 
    and runs a linear mixed effects model. 
    
    Parameters:
        antibody_ids (list): antibody IDs used to find cell types
        pairs_list (list): list containing tuples of each antibody with a cell type
            Example:
                [(antibody_id, idCL), (antibody_id, idCL) ... ]
    Output:
        ab_idCL_lmm_dicts (dict): nested dictionary containing LMM results for
            each cell type for an antibody
            
            Example: 
                {
                  antibody_id_1: {
                      idCL_1: LMM_result_1
                      idCL_2: LMM_result_2
                  },

                  antibody_id_2: {
                      idCL_1: LMM_result_1
                      idCL_2: LMM_result_2
                  },
                }
            
            Where LMM_result is a tuple (LMM Summary table, proportion_background)
    """
    # Create a Spark session
    spark = SparkSession.builder.appName("ParallelMySQLQuery").getOrCreate()

    # We're already provided the idCLs to iterate over
    num_parallel_queries = 100

    # Create an RDD 
    rdd_idCLs = spark.sparkContext.parallelize(pairs_list, num_parallel_queries)

    # Run the spark_query_ab_lmm in parallel and store in a list of dataframes
    results = rdd_idCLs.map(lambda x: spark_query_cell_type_lmm(x[0], x[1])).collect()

    # Stop spark
    spark.stop()

    # Organize results based on antibody
    ab_idCL_lmm_dicts = {}
    for antibody in antibody_ids:
        # Find all entries in results that have this antibody
        # We want the idCL, LMM result and proportion value stored as a tuple
        lmm_res_for_ab = [(x[1], x[2], x[3]) for x in results if x[0] == antibody]

        # Make a temp dict of idCLs
        if lmm_res_for_ab:
            idCLs_dict_temp = {}
            for lmm_res in lmm_res_for_ab:
                idCLs_dict_temp[lmm_res[0]] = (lmm_res[1], lmm_res[2])
        
            ab_idCL_lmm_dicts[antibody] = idCLs_dict_temp
        else:
            ab_idCL_lmm_dicts[antibody] = '{}'

    return ab_idCL_lmm_dicts

def lmm_cell_type_dict_to_df(lmm_results: dict) -> dict:
    """
    Converts the output of the second Spark job in the get_celltypes pipeline
    into a dictionary containing a formatted pandas DataFrame for each cell type

    Parameters:
        lmm_results (dict): dictionary containing LMM results for each cell type
            used for a particular antibody

            Example: 
                {
                  antibody_id_1: {
                      idCL_1: LMM_result_1
                      idCL_2: LMM_result_2
                  },

                  antibody_id_2: {
                      idCL_1: LMM_result_1
                      idCL_2: LMM_result_2
                  },
                }

            Where LMM_result is a tuple (LMM Summary table, proportion_background)
    
    Returns:
        final_ab_lmm_df_dict (dict): dictionary of antibodies with formatted
            dataframes containing LMM results of cell types

            Example: 
                {
                  antibody_id_1: {LMM output (pd.DataFrame)},
                  antibody_id_2: {LMM output (pd.DataFrame)}
                }
    """
    # Dictionary to store all LMM tables, where
    # key: antibody_id, value: LMM table (idCLs are the index)
    final_ab_lmm_df_dict = {}
    
    # Loop over each antibody_id (key) in the outer dictionary
    for antibody_id, idCL_lmm_result in lmm_results.items():
        # EDGE CASE: if there were no lmm results for an antibody, it's value is a STRING called "{}"
        # Append an empty dictionary string to the final dict
        if idCL_lmm_result == "{}":
            final_ab_lmm_df_dict[antibody_id] = "{}"  # String, since it has to be converted to JSON later
            continue
    
        # idCL_lmm_result is another dictionary where keys are cell types (idCL)
        # the value is a tuple(df, float): lmm_output_table, and proportion of background cells
        
        # Start making the output dataframe containing idCLs, for this current antibody_id
        # Index of this dataframe
        idCLs = list(idCL_lmm_result.keys())

        cellType_coeff = []
        stdError = []
        p_values = []
        idCL_dict = {}
        
        # Last row will be for background cell percentage with that idCL and antibody
        idCL_dict["expressed"] = []

        for idCL, table_and_proportion_tuple in idCL_lmm_result.items():
            df = table_and_proportion_tuple[0]
            proportion = table_and_proportion_tuple[1]
            
            # Taking the negative of the coefficient
            cellType_coeff.append(-1 * float(df.loc['cellType']['Coef.']))
            stdError.append(float(df.loc['cellType']['Std.Err.']))
            p_values.append(float(df.loc['cellType']['P>|z|']))

            # Take (1 - proportion_background), which will instead show proportion of signal cells
            idCL_dict["expressed"].append(1 - proportion)

        # Apply Benjamini Hochberg correction to all the p_values
        q_values = multipletests(p_values, method="fdr_bh")[1]

        table_dict = {'coeff': cellType_coeff,
                'stderr': stdError, 
                'p_val': p_values,
                'q_val': q_values}

        entire_dict = {**table_dict, **idCL_dict}
    
        final_df = pd.DataFrame(index=idCLs, data=entire_dict)

        # Add a readable column at the end by doing a query to the database
        params = _config()
        conn = mysql.connector.connect(**params)
        cursor = conn.cursor()
        
        idCL_query = """SELECT cell_types.label
                        FROM cell_types
                        WHERE cell_types.idCL = (%s);"""
        
        # Add cellType column of names from id
        idCL_conv = []
        for idCL in idCLs:
            cursor.execute(idCL_query, (idCL, ))
            idCL_label = cursor.fetchone()[0]
            idCL_conv.append(idCL_label)
        
        # Add column to final df
        final_df.insert(0, 'cellType', idCL_conv)
    
        # Sort by biggest coefficient first
        final_df.sort_values('coeff', ascending=False, inplace=True)

        # Add this dataframe matched with its antibody_id to the final dict
        final_ab_lmm_df_dict[antibody_id] = final_df

    return final_ab_lmm_df_dict

def get_celltypes(ab_ids: list,
                  idBTO: list = None,
                  idExperiment: list = None) -> dict:
    """
    Wrapper function to call the required queries and Spark jobs to find cell types
    that are marked by a provided list of antibodies, optionally restricted by tissue
    type and experiment ID

    Parameters:
        ab_ids (list): antibody IDs used to find cell types
        idBTO (list): tissue type IDs to restrict cell type search results
        idExperiment: experiment IDs to restrict cell type search results
    
    Returns:
        final_result (dict): dictionary of antibodies with formatted
            dataframes containing LMM results of cell types

            Example: 
                {
                  antibody_id_1: {LMM output (pd.DataFrame)},
                  antibody_id_2: {LMM output (pd.DataFrame)}
                }
    """

    # Call FIRST Spark job: for each antibody, find all of its associated idCLs
    ab_idCL_pair_list = run_spark_cell_type_from_ab(ab_ids, idBTO, idExperiment)

    # Call SECOND Spark job: run a LMM for each antibody and its associated idCLs
    ab_idCL_lmm_dicts = run_spark_cell_type_lmm(ab_ids, ab_idCL_pair_list) 

    # Convert results into a dictionary where key: antibodies, value: dataframe
    final_result = lmm_cell_type_dict_to_df(ab_idCL_lmm_dicts)

    if len(final_result) == 0:
        return ("No cell types found. Please try with different parameters.")

    return final_result

def plot_celltypes(ab_id: str, 
                   id_CLs: list) -> pd.DataFrame:
    """
    Queries a database to find normalized values based on cells found in
    an antibody and certain cell types

    Processes queried data into descriptive statistics used for plotting
    on the client.

    Parameters:
        ab_id (str): antibody ID used for querying
        id_CLs (list): list of cell type IDs used for querying

    Returns:
        summary_stats (pd.DataFrame): dataframe containing the min, max, 
            mean, etc for a particular cell type's normalized values
    """
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)
    
    idCL_placeholders = ','.join(['%s'] * len(id_CLs))
    
    with_idCL_query = """SELECT cells.idCell, cells.idCellOriginal, 
                            antigen_expression.normValue, cells.idExperiment, antigen_expression.idAntibody, cells.idCL
                        FROM antigen_expression
                        INNER JOIN cells ON antigen_expression.idCell = cells.idCell
                        WHERE antigen_expression.idAntibody = (%s)
                        AND cells.idCL IN (%s);""" % ('%s', idCL_placeholders)

    parameters = [ab_id]
    parameters.extend(id_CLs)
    
    cells_with_idCL_df = pd.read_sql(sql=with_idCL_query, params=parameters, con=conn)

    grouped_normvalue_by_celltype = group_normvalue("idCL", cells_with_idCL_df)
    summary_stats = summary_statistics("idCL", grouped_normvalue_by_celltype)

    if conn is not None:
        conn.close()
        
    return summary_stats

#------------------------ Functions for get_experiment ------------------------#
def get_experiments(ab_id: list, 
                    idCL: list = None,
                    idBTO: list = None) -> pd.DataFrame:
    """
    Queries a database to find experiments used based on antibody,
    cell type, or tissue type (includes combinations of parameters)

    Parameters:
        ab_id (list): finds experiments with these antibody IDs
        idCL (list): finds experiments with these cell type IDs
        idBTO (list): finds experiments with these tissue type IDs
    
    Returns:
        exp_df (pd.DataFrame): dataframe containing list of experiments
            returned from the query
    """

    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)
    
    # ab_id
    if ab_id is not None and idCL is None and idBTO is None:
        ab_placeholders = ','.join(['%s'] * len(ab_id))
        ab_query = """SELECT *
                    FROM experiments
                    WHERE experiments.idExperiment IN 
                        (SELECT DISTINCT antigen_expression.idExperiment
                         FROM antigen_expression
                         WHERE antigen_expression.idAntibody IN (%s));""" % (ab_placeholders)
        
        parameters = ab_id.copy()
        exp_df = pd.read_sql(sql=ab_query, params=parameters, con=conn)

    # ab_id and idCL
    elif ab_id is not None and idCL is not None and idBTO is None:
        ab_placeholders = ','.join(['%s'] * len(ab_id))
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
        
    # ab_id and idBTO
    elif ab_id is not None and idCL is None and idBTO is not None:
        ab_placeholders = ','.join(['%s'] * len(ab_id))
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

    # ab_id and idCL and idBTO
    elif ab_id is not None and idCL is not None and idBTO is not None:
        ab_placeholders = ','.join(['%s'] * len(ab_id))
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
    
    # idCL
    elif ab_id is None and idCL is not None and idBTO is None:
        idCLs_placeholders = ','.join(['%s'] * len(idCL))
        idCL_query = """SELECT * 
                        FROM experiments
                        WHERE experiments.idExperiment IN 
                            (SELECT DISTINCT cells.idExperiment
                            FROM cells
                            WHERE cells.idCL IN (%s) );""" % (idCLs_placeholders)
        parameters = idCL.copy()
        exp_df = pd.read_sql(sql=idCL_query, params=parameters, con=conn)
    
    # idCL and idBTO
    elif ab_id is None and idCL is not None and idBTO is not None:
        idCLs_placeholders = ','.join(['%s'] * len(idCL))
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
        
        no_ab_idCL_BTO_query = """SELECT * 
                                FROM experiments
                                WHERE experiments.idExperiment IN 
                                    (SELECT DISTINCT cells.idExperiment
                                    FROM cells
                                    WHERE cells.idCL IN (%s) 
                                    AND experiments.idBTO IN (%s) );""" % (idCLs_placeholders, idBTO_placeholders)
        parameters = idCL.copy()
        parameters.extend(idBTO)
        exp_df = pd.read_sql(sql=no_ab_idCL_BTO_query, params = parameters, con=conn)
    
    # idBTO
    elif ab_id is None and idCL is None and idBTO is not None:
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
        idBTO_query = """SELECT * 
                        FROM experiments
                        WHERE experiments.idBTO IN (%s);""" % (idBTO_placeholders)
        
        parameters = idBTO.copy()
        exp_df = pd.read_sql(sql=idBTO_query, params=parameters, con=conn)
    
    if len(exp_df.index) == 0:
        return ("No experiments matching that search query were found in the database.")

    # Add tissue name of idBTO
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()
    
    tissue_query = """SELECT tissues.tissueName
                    FROM tissues
                    WHERE tissues.idBTO = (%s);"""
    
    # Add cellType column of names from id
    tissue_conv = []
    for idBTO in exp_df['idBTO']:
        cursor.execute(tissue_query, (idBTO, ))
        idCL_label = cursor.fetchone()[0]
        tissue_conv.append(idCL_label)
    
    # Add column to final df
    exp_df['tissue'] = tissue_conv

    if conn is not None:
        cursor.close()
        conn.close()

    return exp_df

#------------------------ Functions for which_antibodies ------------------------#
def get_unique_antibodies() -> list:
    """
    Finds all unique antibodies in the database

    Parameters:
        N/A
    
    Returns:
        ab_list (list): antibodies currently in the database
    """
    warnings.filterwarnings('ignore')

    
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()
    
    ab_query = """SELECT antibodies.idAntibody
                 FROM antibodies;"""
    
    cursor.execute(ab_query)
    res = cursor.fetchall()
    ab_list = [ab[0] for ab in res]
    
    if conn is not None:
        cursor.close()
        conn.close()
    
    return ab_list

def which_antibodies(search_query: str) -> pd.DataFrame:
    """
    Performs a full-text search query in a database to find
    the most relevant antibodies based on user input
    
    Parameters:
        search_query (str): a phrase or word containing the desired 
            antibody
    
    Returns:
        final_df (pd.DataFrame): dataframe containing best matching 
            results for the provided query 
    """
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)

    which_ab_query = """SELECT * 
                    FROM antibodies
                    WHERE MATCH (idAntibody, host, cloneID, abName, abTarget, citation, comments, vendor, catalogNum)
                    AGAINST ( (%s) IN NATURAL LANGUAGE MODE);""" % ('%s')

    parameters = [search_query]
    exp_df = pd.read_sql(sql=which_ab_query, params=parameters, con=conn)

    # If none of the antibodies from the query were found in the database
    if len(exp_df.index) == 0:
        return ("No antibodies matching that search query were found in the database.")

    # Reorder columns
    exp_df = exp_df[['idAntibody', 'abName', 'abTarget', 
                     'polyclonal', 'citation', 'comments', 
                     'cloneID', 'host', 'vendor', 'catalogNum']]
    
    # Find and replace polyclonal column with appropriate clonality
    exp_df.loc[exp_df['polyclonal'] == 0, 'polyclonal'] = 'monoclonal'
    exp_df.loc[exp_df['polyclonal'] == 1, 'polyclonal'] = 'polyclonal'

    final_df = exp_df.rename(columns={'polyclonal': 'clonality'})

    experiment_each_ab = []
    for ab in list(final_df['idAntibody']):
        
        res_df = get_experiments(ab_id=[ab])
        # Edge case: if there are no experiments for a result, it will return string
        if isinstance(res_df, str):
            # Append None in the result column
            experiment_each_ab.append("None")
        else:
            experiments_found = list(res_df['idExperiment'])

            # Convert all experiments found into a string
            exp_found_str = ','.join(str(exp) for exp in experiments_found)
            experiment_each_ab.append(exp_found_str)

    # Add to final column
    final_df['idExperiment_used'] = experiment_each_ab

    if conn is not None:
        conn.close()

    return final_df


#------------------------ Functions for which_celltypes ------------------------#
def get_unique_idCLs() -> list:
    """
    Find all unique cell types in the database

    Parameters:
        N/A

    Returns:
        idCLs_list (list): cell types in the database
    """
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()
    
    idCLs_query = """SELECT cell_types.idCL 
                    FROM cell_types;"""
    
    cursor.execute(idCLs_query)
    res = cursor.fetchall()
    idCLs_list = [idCL[0] for idCL in res]
    
    if conn is not None:
        cursor.close()
        conn.close()
    
    return idCLs_list

def get_idCL_experiments(idCL: str) -> pd.DataFrame:
    """
    Finds all experiments used for a provided cell type

    Parameters:
        idCL (str): cell type ID
    
    Returns:
        idCL_exp_df (pd.DataFrame): dataframe containing list of experiments
            for a cell type
    """
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()
    
    idCL_exp_query = """SELECT *
                    FROM experiments
                    WHERE experiments.idExperiment IN
                        (SELECT DISTINCT cells.idExperiment 
                        FROM cells
                        WHERE cells.idCL IN ('%s') );""" % (idCL)
    
    idCL_exp_df = pd.read_sql(sql=idCL_exp_query, con=conn)
    
    if conn is not None:
        cursor.close()
        conn.close()
    
    return idCL_exp_df

def which_celltypes(search_query:str) -> pd.DataFrame:
    """
    Performs a full-text search query in a database to find
    the most relevant cell types based on user input
    
    Parameters:
        search_query (str): a phrase or word containing the desired 
            cell type
    
    Returns:
        final_df (pd.DataFrame): dataframe containing best matching 
            results for the provided query 
    """  
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)
    
    ct_query = """SELECT cell_types.idCL, cell_types.label  
                FROM cell_types
                WHERE MATCH (idCL, label)
                AGAINST ( (%s) IN NATURAL LANGUAGE MODE);""" % ('%s')
    
    parameters = [search_query]
    ct_exp_df = pd.read_sql(sql=ct_query, params=parameters, con=conn) 
    
    # FIRST CHECK: Are there any results from the full text search query alone?
    if len(ct_exp_df.index) == 0:
        return ("No cell types matching that search query were found in the database.")

    # Cleaned input:
    cleaned_input = search_query.lower().strip()
    split_cleaned = re.split("\W+", cleaned_input)

    final_df = ct_exp_df[ct_exp_df['label'].str.lower()
                         .apply(lambda sentence: all(word in re.split("\W+", sentence) for word in split_cleaned))]
    
    # SECOND CHECK: Are there any results after we filter the results?
    if len(final_df.index) == 0:
        return ("No cell types matching that search query were found in the database.")

    experiment_each_idCL = []
    # Find experiments using the idCL
    for idCL in list(final_df['idCL']):
        res_df = get_idCL_experiments(idCL=idCL)
        experiments_found = list(res_df['idExperiment'])

        # Convert all experiments found into a string
        exp_found_str = ','.join(str(exp) for exp in experiments_found)
        experiment_each_idCL.append(exp_found_str)
    
    # Add final column indicating the experiments used
    final_df['idExperiment_used'] = experiment_each_idCL

    # Reset the index after we're done filtering
    final_df.reset_index(drop=True, inplace=True)
    
    if conn is not None:
        conn.close()

    return final_df

#------------------------ Functions for which_experiments ------------------------#
def which_experiments(search_query: str):
    """
    Performs a partial-match search query in a database to find
    the most relevant experiments based on user input
    
    Parameters:
        search_query (str): a phrase or word containing the desired 
            experiment
    
    Returns:
        final_df (pd.DataFrame): dataframe containing best matching 
            results for the provided query 
    """
    warnings.filterwarnings('ignore')
    
    params = _config()
    conn = mysql.connector.connect(**params)

    exp_query = """SELECT e.idExperiment, e.nameExp, e.typeExp, e.pmid, e.doi, e.idBTO, t.tissueName
            FROM experiments e
            INNER JOIN tissues t ON e.idBTO = t.idBTO
            WHERE LOWER(CONCAT_WS(' ', e.idExperiment, e.nameExp, e.typeExp, e.pmid, e.doi, e.idBTO, t.tissueName)) LIKE CONCAT('%%', %s, '%%')
            ORDER BY e.idExperiment ASC"""
    
    parameters = [search_query]
    exp_df = pd.read_sql(sql=exp_query, params=parameters, con=conn)
    exp_df.rename(columns={"tissueName": "tissue"}, inplace=True)

    if len(exp_df.index) == 0:
        return ("No experiments matching that search query were found in the database.")

    return exp_df

#------------------------ Functions for querying reference data for STvEA ------------------------#
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

def _uniprot_aliases(sci_crunch_alias: str = None, 
                     user_uniprotID: str = None) -> tuple:
    """
    Retrieves information about an antibody using UniProt's API

    Parameters:
        sci_crunch_alias (str): antibody alias returned from SciCrunch API
        uniprotID (str): user-provided UniProtID for manual entry

    Returns:
        uniprotID, otherAliases (tuple [str, list]): antibody alias
            and aliases retrieved from UniProt, or an empty list for no results
    """

    otherAliases = []

    # Removed homo sapiens organism id filter from all api queries
                
    # Protein and gene search (strictest level)
    protein_gene_params = {
        'query':  f'protein_name:{sci_crunch_alias} AND gene_exact:{sci_crunch_alias} AND reviewed:true',
        'fields': 'gene_primary, gene_synonym, protein_name',
        'format': 'json'
    }

    # Protein search query params
    protein_params = {
        'query':  f'protein_name:{sci_crunch_alias} AND reviewed:true',
        'fields': 'gene_primary, gene_synonym, protein_name',
        'format': 'json'
    }

    # Gene search query params if protein query fails (response != 200 or no results)
    gene_params = {
        'query':  f'gene_exact:{sci_crunch_alias} AND reviewed:true',
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
        # print("resJSON from uniprot pair matching:", resJSON)
        uniprotID = resJSON['results'][0]['primaryAccession']
    except:
        return []
        # raise Exception(f"Unable to find UniProtID for: {sci_crunch_alias}")

    # If there was a UniprotID, we can get its aliases
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
    
def get_cloneID(ab_id: str) -> str:
    """
    Retrieves cloneID for an antibody ID from SciCrunch
    """
    ### SciCrunch ###
    sciJSON = _sci_crunch_hits(ab_id)

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
            return cloneID
    
def ab_subset(idBTO: list = None, idExperiment: list = None) -> list:
    """ 
    Used for narrowing down the search space for antibodies. If this is called, then
    we can eliminate some of the antibodies from the spreadsheet to save redundant queries

    Then, we do further analysis for all the antibodies not shared
    """
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()

    similar_abs_list = []

    # Filter by idBTO
    if idBTO is not None and idExperiment is None:
        # Generate placeholders
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
        
        ab_idBTO_query = """SELECT DISTINCT(antibodies.idAntibody)
                    FROM cells
                    INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
                    INNER JOIN tissues on experiments.idBTO = tissues.idBTO
                    INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                    INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
                    WHERE experiments.idBTO IN (%s);""" % (idBTO_placeholders)
        
        parameters = idBTO
        cursor.execute(ab_idBTO_query, parameters)
        similar_abs = cursor.fetchall()
        similar_abs_list = [ab[0] for ab in similar_abs]

    # Filter by idExperiment
    elif idBTO is None and idExperiment is not None:
        idExp_placeholders = ','.join(['%s'] * len(idExperiment))
        
        all_idExp_query = """SELECT DISTINCT(antibodies.idAntibody)
                    FROM cells
                    INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
                    INNER JOIN tissues on experiments.idBTO = tissues.idBTO
                    INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                    INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
                    WHERE experiments.idExperiment IN (%s);""" % (idExp_placeholders)
        
        parameters = idExperiment
        cursor.execute(all_idExp_query, parameters)
        similar_abs = cursor.fetchall()
        similar_abs_list = [ab[0] for ab in similar_abs]

    elif idBTO is not None and idExperiment is not None:
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
        idExp_placeholders = ','.join(['%s'] * len(idExperiment))
        
        all_idBTO_idExp_query = """SELECT DISTINCT(antibodies.idAntibody)
                            FROM cells
                            INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
                            INNER JOIN tissues on experiments.idBTO = tissues.idBTO
                            INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                            INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
                            WHERE experiments.idBTO IN (%s) AND
                            	  experiments.idExperiment IN (%s);""" % (idBTO_placeholders, idExp_placeholders)
        
        parameters = idBTO
        parameters.extend(idExperiment)
    
        cursor.execute(all_idBTO_idExp_query, parameters)
        similar_abs = cursor.fetchall()
        similar_abs_list = [ab[0] for ab in similar_abs]

    return similar_abs_list

def extract_alnum_lowercase(string) -> str:
    """
    Extracts all numbers and letters from a string as lowercase

    Parameters:
        string (str): character input
    
    Returns:
        string (str): filtered characters from string
    """
    return ''.join(ch for ch in string if ch.isalnum()).lower()

def match_antibody(ab_name: str, 
                   ab_id: str, 
                   option: int = 1, 
                   idBTO: list = None, 
                   idExperiment: list = None) -> list:
    """
    Finds the most relevant antibodies in the database given an antibody name and ID.
    If an exact antibody is not found in the database, a matching is performed that
    considers the clone ID and alias for that antibody to find the next best result.
    Allows different levels of search/strictness.

    Parameters:
        ab_name (str): antibody name provided on the experiment spreadsheet
        ab_id (str): antibody ID provided on the experiment spreadsheet
        option (int): level of strictness when searching, where
            1: parse by clone ID and alias (default)
            2: parse by alias and antibody ID (most relaxed)
            3: parse by antibody ID (strictest)
    
    Returns:
        matched_antibodies (list): list of antibodies returned from the database
            matching the provided input antibody from the experiment spreadsheet

    """    
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()

    # If the antibody isn't in the database, start trying to match the next best result
    # Default option: match based on cloneID and aliases
    if option == 1:
        aliases = []

        # Get clone ID 
        cloneID = get_cloneID(ab_id)

        # Get aliases if there is also a valid UniProtID for this antibody 
        # If there are none, "aliases" remains as an empty list []
        try:
            uniprot, found_aliases = _uniprot_aliases(sci_crunch_alias=ab_name)
            aliases = found_aliases
        except:
            aliases = []

        # Holds all similar antibodies
        similar_abs_list = []
        
        # Find all rows (idAntibody) in the database that share the same clone ID (with varying filters)
        # For only clone ID
        if idBTO is None and idExperiment is None:
            find_ab_similar_clone_query = """SELECT antibodies.idAntibody 
                                          FROM antibodies
                                          WHERE antibodies.cloneID = (%s);"""
    
            cursor.execute(find_ab_similar_clone_query, (cloneID, ))
            similar_abs = cursor.fetchall()
            similar_abs_list = [ab[0] for ab in similar_abs]

        # For clone ID and tissue filter
        elif idBTO is not None and idExperiment is None:
            find_ab_tissue = """SELECT DISTINCT antibodies.idAntibody
                        FROM cells
                        INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
                        INNER JOIN tissues ON experiments.idBTO = tissues.idBTO
                        INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                        INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
                        WHERE antibodies.cloneID = %s AND
                              experiments.idBTO IN (%s);"""

            parameters = [cloneID]
            parameters.extend(idBTO)
        
            cursor.execute(find_ab_tissue, parameters)
            similar_abs = cursor.fetchall()
            similar_abs_list = [ab[0] for ab in similar_abs]

        # For clone ID and idExperiment filter
        elif idBTO is None and idExperiment is not None:
            find_ab_exp = """SELECT DISTINCT (antibodies.idAntibody)
                        FROM cells
                        INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
                        INNER JOIN tissues on experiments.idBTO = tissues.idBTO
                        INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                        INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
                        WHERE antibodies.cloneID = %s AND
                        experiments.idExperiment IN (%s);"""

            parameters = [cloneID]
            parameters.extend(idExperiment)
        
            cursor.execute(find_ab_exp, parameters)
            similar_abs = cursor.fetchall()
            similar_abs_list = [ab[0] for ab in similar_abs]

        # For clone ID and tissue AND idExperiment filter
        elif idBTO is not None and idExperiment is not None:
            find_ab_tissue_exp = """SELECT DISTINCT (antibodies.idAntibody)
                                FROM cells
                                INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
                                INNER JOIN tissues on experiments.idBTO = tissues.idBTO
                                INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                                INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
                                WHERE antibodies.cloneID = %s AND
                                experiments.idBTO IN (%s) AND
                                experiments.idExperiment IN (%s);"""

            parameters = [cloneID]
            parameters.extend(idBTO)
            parameters.extend(idExperiment)
        
            cursor.execute(find_ab_tissue_exp, parameters)
            similar_abs = cursor.fetchall()
            similar_abs_list = [ab[0] for ab in similar_abs]

        # If there were more than 1 'similar' antibody based on cloneID
        # Investigate further by checking for a matching alias between 
        # antibody in database and the unknown antibody. 
        # This will only happen if there was a valid UniProt ID found earlier
        if len(similar_abs_list) > 0 and len(aliases) > 0:

            matched_antibodies = []
            
            for similar_ab in similar_abs_list:
                # Find aliases for the antibody in the database
                aliases_query = """SELECT aliases.alias FROM antigens
                                INNER JOIN aliases on aliases.idUniProtKB = antigens.idUniProtKB
                                WHERE antigens.idAntibody = (%s);"""
                cursor.execute(aliases_query, (similar_ab, ))
                similar_aliases = cursor.fetchall()
                similar_aliases_list = [alias[0] for alias in similar_aliases]

                # Find any similar results
                # Remove any symbols and convert all strings to lowercase
                aliases_modified = [extract_alnum_lowercase(x) for x in aliases]
                similar_aliases_list_modified = [extract_alnum_lowercase(x) for x in similar_aliases_list]
                common_aliases = [i for i in aliases_modified if i in similar_aliases_list_modified]

                if len(common_aliases) > 0:
                    matched_antibodies.append(similar_ab)
                else:
                    continue # skip to next result

            return matched_antibodies

        # If only searching based on cloneID
        elif len(similar_abs_list) > 0:
            return similar_abs_list

        # No results from matching by cloneID and/or aliases
        else:
            return []

    # Option 2: Search by only alias and antibody target, no clone ID
    # If there are no aliases (no valid UniprotID), then only search by antibody target
    # if ab_exists_result != 1 and option == 2:
    if option == 2:
        aliases = []

        # Get aliases if there is also a valid UniProtID for this antibody 
        # If there are none, "aliases" remains as an empty list []
        try:
            uniprot, found_aliases = _uniprot_aliases(sci_crunch_alias=ab_name)
            aliases = found_aliases
        except:
            aliases = []

        # Holds all similar antibodies
        similar_abs_list = []

        # Find all rows (idAntibody) in the database that share the same abTarget/ab_name (with varying filters)
        # For only ab_name
        if idBTO is None and idExperiment is None:
            # Find all rows (idAntibody) in the database that share the same abTarget
            find_ab_similar_target_query = """SELECT antibodies.idAntibody 
                                          FROM antibodies
                                          WHERE antibodies.abTarget = (%s);"""
    
            cursor.execute(find_ab_similar_target_query, (ab_name, ))
            similar_abs = cursor.fetchall()
            similar_abs_list = [ab[0] for ab in similar_abs]

        # For ab_name and tissue filter
        elif idBTO is not None and idExperiment is None:
            find_ab_tissue = """SELECT DISTINCT antibodies.idAntibody
                        FROM cells
                        INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
                        INNER JOIN tissues ON experiments.idBTO = tissues.idBTO
                        INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                        INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
                        WHERE antibodies.abTarget = %s AND
                              experiments.idBTO IN (%s);"""

            parameters = [ab_name]
            parameters.extend(idBTO)
        
            cursor.execute(find_ab_tissue, parameters)
            similar_abs = cursor.fetchall()
            similar_abs_list = [ab[0] for ab in similar_abs]

        # For ab_name and idExperiment filter
        elif idBTO is None and idExperiment is not None:
            find_ab_exp = """SELECT DISTINCT antibodies.idAntibody
                        FROM cells
                        INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
                        INNER JOIN tissues ON experiments.idBTO = tissues.idBTO
                        INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                        INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
                        WHERE antibodies.abTarget = %s AND
                        experiments.idExperiment IN (%s);"""

            parameters = [ab_name]
            parameters.extend(idExperiment)
        
            cursor.execute(find_ab_exp, parameters)
            similar_abs = cursor.fetchall()
            similar_abs_list = [ab[0] for ab in similar_abs]

        # For ab_name and tissue AND idExperiment filter
        elif idBTO is not None and idExperiment is not None:
            find_ab_tissue_exp = """SELECT DISTINCT (antibodies.idAntibody)
                                FROM cells
                                INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
                                INNER JOIN tissues on experiments.idBTO = tissues.idBTO
                                INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
                                INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
                                WHERE antibodies.abTarget = %s AND
                                experiments.idBTO IN (%s) AND
                                experiments.idExperiment IN (%s);"""

            parameters = [ab_name]
            parameters.extend(idBTO)
            parameters.extend(idExperiment)
        
            cursor.execute(find_ab_tissue_exp, parameters)
            similar_abs = cursor.fetchall()
            similar_abs_list = [ab[0] for ab in similar_abs]

        # If there were more than 1 'similar' antibody based on cloneID
        # Investigate further by checking for a matching alias between 
        # antibody in database and the unknown antibody
        if len(similar_abs_list) > 0 and len(aliases) > 0:
            matched_antibodies = []
            
            for similar_ab in similar_abs_list:
                # Find aliases for the antibody in the database
                aliases_query = """SELECT aliases.alias FROM antigens
                                INNER JOIN aliases on aliases.idUniProtKB = antigens.idUniProtKB
                                WHERE antigens.idAntibody = (%s);"""
                cursor.execute(aliases_query, (similar_ab, ))
                similar_aliases = cursor.fetchall()
                similar_aliases_list = [alias[0] for alias in similar_aliases]

                # Find any similar results
                # Remove any symbols and convert all strings to lowercase
                aliases_modified = [extract_alnum_lowercase(x) for x in aliases]
                similar_aliases_list_modified = [extract_alnum_lowercase(x) for x in similar_aliases_list]
                common_aliases = [i for i in aliases_modified if i in similar_aliases_list_modified]

                if len(common_aliases) > 0:
                    matched_antibodies.append(similar_ab)
                else:
                    continue # skip to next result

            return matched_antibodies
            
        # If only searching based on abTarget
        elif len(similar_abs_list) > 0:
            return similar_abs_list

        # If no results from abTarget and aliases
        else:
            return []

    # Option 3: strict search by antibody ID only in the database
    if option == 3:
        check_ab_exists_query = """SELECT antibodies.idAntibody
                                FROM antibodies 
                                WHERE antibodies.idAntibody=(%s);"""

        cursor.execute(check_ab_exists_query, (ab_id, ))
        ab_exists_result = cursor.fetchone()
        
        if ab_exists_result:
            return [ab_id]
        else:
            return []
        
def parse_antibodies(antibody_pairs: list, option = 1, idBTO=None, idExperiment=None) -> dict:
    """
    Parses each antibody name and ID in the experiment spreadsheet to 
    find matching antibodies in the database.

    Parameters:
        antibody_pairs (list): list containing: [antibody_name, antibody_id]
        option (int): degree of strictness for matching
    
    Returns:
        database_to_original_ab_dict (dict): dictionary containing IDs where:
            key: antibody ID from database, value: original antibody ID
    """
    antibodies_to_query = []
    database_to_original_ab_dict = {} # keys: ab retrieved from db, value: original antibody from dataset

    # Get all antibodies ID from the spreadsheet
    antibody_in_spreadsheet = [antibody_and_id[1] for antibody_and_id in antibody_pairs]

    # Narrow the ab search space if provided idBTO or idExp
    if idBTO is not None or idExperiment is not None:
        # Find all antibodies from the idBTO/idExperiment
        ab_in_db = ab_subset(idBTO=idBTO, idExperiment=idExperiment)
    
        # Find shared antibodies
        shared_antibodies = list(set(antibody_in_spreadsheet) & set(ab_in_db))
    
        # For all shared antibodies, add them to the dict. Remove this from ab in spreadsheet to reduce the number of queries
        for ab in shared_antibodies:
            database_to_original_ab_dict[ab] = ab
            antibody_in_spreadsheet.remove(ab)

    # At this point, we can start performing the queries for the remaining antibodies
    remaining_ab_pairs = [antibody_and_id for antibody_and_id in antibody_pairs if antibody_and_id[1] in antibody_in_spreadsheet]

    for antibody_and_id in remaining_ab_pairs:
        print(f"Looking up: {antibody_and_id[0]}, {antibody_and_id[1]}")
        result = match_antibody(antibody_and_id[0], antibody_and_id[1], option=option, idBTO=idBTO, idExperiment=idExperiment)
        if len(result) > 0: # If multiple antibodies due to matching, take first hit
            antibodies_to_query.append(result[0])
            database_to_original_ab_dict[result[0]] = antibody_and_id[1]
        else:
            continue

    return database_to_original_ab_dict

def drop_antibodies(truth_table: pd.DataFrame) -> pd.DataFrame:
    """
    Removes columns (antibodies) from a table that contains all 0s. This table contains
    boolean integers (0 = false, 1 = true) where rows = experiment ID and columns = antibody ID.
    If an antibody is not used in all experiments, it is dropped/removed from the table.

    Parameters:
        truth_table (pd.DataFrame): table containing a 0 or 1 depending on expression
            of a given antibody in an experiment
        
    Returns:
        truth_final (pd.DataFrame): table containing only antibodies that are used
            in at least one experiment
    """
    columns_to_exclude = ['idExperiment']
    temp_columns = truth_table["idExperiment"]
    temp_truth = truth_table.loc[:, ~truth_table.columns.isin(columns_to_exclude)]

    # Remove any multi index
    truth_table.columns.name = None

    # Perform first antibody drop here, for all columns with entire values set to 0
    # Check if all values in each column are 0
    zero_columns = (temp_truth == 0).all()
    
    # Get the column names where all values are 0
    columns_to_drop = zero_columns[zero_columns].index
    columns_to_drop
    
    # # Drop the columns with all 0s
    truth_final = temp_truth.drop(columns=columns_to_drop)

    # Restore the rest of the non-antibody columns
    truth_final.insert(0, 'idExperiment', temp_columns)
    
    return truth_final

def find_exp_truth_table(antibodies: list,
                         idBTO: list = None,
                         idExperiment: list = None) -> pd.DataFrame:
    """
    Given a list of antibodies, find where each antibody is used in all the
    experiments in the database. Optionally filter by a specific tissue type
    or particular experiments.

    Parameters:
        antibodies (list): list of antibody IDs
        idBTO (list): list of tissue IDs to filter (optional)
        idExperiment (list): list of experiment IDs to filter (optional)
    
    Returns:    
        exp_table (pd.DataFrame): table where rows = experiment IDs and 
            columns = antibody IDs. Value of each cell can be 0 or 1, where
            0 indicates that the antibody is not used in that experiment and
            1 indicates that the antibody is used in the experiment
    """
    
    warnings.filterwarnings('ignore')
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()

    # No additional filtering
    if antibodies is not None and idBTO is None and idExperiment is None:
        ab_placeholders = ','.join(['%s'] * len(antibodies))
                       
        ab_query = """SELECT DISTINCT (antibodies.idAntibody), antibodies.abTarget, experiments.idExperiment
            FROM cells
            INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
            INNER JOIN tissues on experiments.idBTO = tissues.idBTO
            INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
            INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
            WHERE antibodies.idAntibody IN (%s);""" % (ab_placeholders)
    
        parameters = antibodies.copy()
        
        # Raw data from database 
        # where each row: [idAntibody, abTarget, idExperiment]
        
        exp_table = pd.read_sql(sql=ab_query, params=parameters, con=conn)
                             
     # Filter only by experiment ID
    elif antibodies is not None and idBTO is None and idExperiment is not None:
        ab_placeholders = ','.join(['%s'] * len(antibodies))
        idExperiment_placeholders = ','.join(['%s'] * len(idExperiment))
                       
        ab_idExp_query = """SELECT DISTINCT (antibodies.idAntibody), antibodies.abTarget, experiments.idExperiment
            FROM cells
            INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
            INNER JOIN tissues on experiments.idBTO = tissues.idBTO
            INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
            INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
            WHERE antibodies.idAntibody IN (%s) AND
             experiments.idExperiment IN (%s);""" % (ab_placeholders,  idExperiment_placeholders)
    
        parameters = antibodies.copy()
        parameters.extend(idExperiment)
    
        exp_table = pd.read_sql(sql=ab_idExp_query, params=parameters, con=conn)

    # Filter only by tissue ID
    elif antibodies is not None and idBTO is not None and idExperiment is None:
        ab_placeholders = ','.join(['%s'] * len(antibodies))
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
                       
        ab_idBTO_query = """SELECT DISTINCT (antibodies.idAntibody), antibodies.abTarget, experiments.idExperiment
            FROM cells
            INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
            INNER JOIN tissues on experiments.idBTO = tissues.idBTO
            INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
            INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
            WHERE antibodies.idAntibody IN (%s) AND
             experiments.idBTO IN (%s);""" % (ab_placeholders, idBTO_placeholders)
    
        parameters = antibodies.copy()
        parameters.extend(idBTO)
        
        exp_table = pd.read_sql(sql=ab_idBTO_query, params=parameters, con=conn)

    # Filter by both tissue and experiment
    elif antibodies is not None and idBTO is not None and idExperiment is not None:
        ab_placeholders = ','.join(['%s'] * len(antibodies))
        idBTO_placeholders = ','.join(['%s'] * len(idBTO))
        idExperiment_placeholders = ','.join(['%s'] * len(idExperiment))
                       
        ab_idBTO_idExp_query = """SELECT DISTINCT (antibodies.idAntibody), antibodies.abTarget, experiments.idExperiment
            FROM cells
            INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
            INNER JOIN tissues on experiments.idBTO = tissues.idBTO
            INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
            INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
            WHERE antibodies.idAntibody IN (%s) AND
             experiments.idBTO IN (%s) AND 
             experiments.idExperiment IN (%s);""" % (ab_placeholders, idBTO_placeholders, idExperiment_placeholders)
    
        parameters = antibodies.copy()
        parameters.extend(idBTO)
        parameters.extend(idExperiment)

        exp_table = pd.read_sql(sql=ab_idBTO_idExp_query, params=parameters, con=conn)
        
    # Format table into rows: experiments, columns: antibodies, cell value: abTarget
    exp_pivot = exp_table.pivot_table(index='idExperiment', columns='idAntibody', 
                                      values='abTarget', aggfunc='first')

    # Create truth table, where 0: not expressed, 1: expressed (IOW: there's an abTarget present)
    truth_table = (~exp_pivot.isna()).astype(int)

    # Remove any multi index
    truth_table.columns.name = None
                             
    # Move idExperiment into a column, and reset the index so it starts counting at 0, 1, 2...
    # Reason being that the experiment IDs will not always be in a sequence
    # Index will have to be matched back up to the value in column 'idExperiment'
    
    pivot_table_reset = truth_table.reset_index()
    pivot_table_reset.index = range(len(pivot_table_reset))         

    # table is ready for pairwise calculations
    return pivot_table_reset

def pairwise_distance_matrix(exp_ab_truth_table: pd.DataFrame) -> pd.DataFrame:
    """
    Given the experiment and antibody truth table (0s and 1s), calculate a 
    pairwise distance metric used to measure the number of antibodies
    associated to a particular experiment 

    Parameters:
        exp_ab_truth_table (pd.DataFrame): experiment and antibody
            truth table for which antibodies are expressed in each experiment

    Returns:
        distance_matrix (pd.DataFrame): pairwise distance matrix for
            every experiment and antibody
    """
    num_exp = len(exp_ab_truth_table.index)
    num_ab = len(exp_ab_truth_table.columns) - 1 # disregard idExperiment
    
    distance_matrix = np.zeros([num_exp, num_exp])

    for i in range(0, num_exp): ## need to start at 1 due to exp ID
        for j in range(i + 1, num_exp):
            ab_a = list(exp_ab_truth_table.loc[i, exp_ab_truth_table.columns != 'idExperiment'])
            ab_b = list(exp_ab_truth_table.loc[j, exp_ab_truth_table.columns != 'idExperiment'])

            ab_bool = [x == 1 and y == 1 for x, y in zip(ab_a, ab_b)]
            ab_pair_bool = list(np.where(ab_bool)[0]) # find number of pairs of exp that have info about an antibody
            num_exp_with_info_about_ab = len(ab_pair_bool)
            
            # Calculate Distance Metric
            metric = (num_ab - num_exp_with_info_about_ab) / (num_ab)
                
            # Add to distance matrix
            distance_matrix[i][j] = metric
            distance_matrix[j][i] = metric
            
    # Matrix is ready to be converted into adjacency matrix
    return distance_matrix

def pairwise_adjacency_matrix(pairwise_distance_matrix: pd.DataFrame,
                              filter: float = 1) -> pd.DataFrame:
    """
    Applies a filter to the pairwise distance matrix and returns a 
    boolean (0 = False, 1 = True) for each value. The resulting boolean
    matrix is used as an adjacency matrix for a graph representing 
    the relationship between each experiment (node). 

    Parameters:
        pairwise_distance_matrix (pd.DataFrame): pairwise distance 
            matrix for every experiment and antibody
        filter (float): threshold to decide whether a pairwise distance
            can be used to account for a given experiment
    
    Returns:
        adjacency_matrix (pd.DataFrame): adjacency matrix consisting of
            0s and 1s after filtering by a threshold.
    """
    # Create a boolean mask for values under a threshold
    # We will only consider pairs that fall under this threshold
    filtered_bool_mask = pd.DataFrame(pairwise_distance_matrix) < filter

    # Create adjacency matrix where 0 = false, 1 = true
    adjacency_matrix = filtered_bool_mask.astype(int)
    return adjacency_matrix

def find_giant_component(adjacency_matrix: pd.DataFrame) -> list:
    """
    Finds the giant component in a graph of nodes representing
    experiment IDs. This component contains the experiments that will be
    considered when querying for the reference dataset. The purpose of this
    giant component is to reduce the chance of querying data from experiments
    that have no expression of an antibody of interest. 

    Parameters:
        adjacency_matrix (pd.DataFrame): matrix containing 0s and 1s
            used to create a graph representing experiments
    
    Returns:
        G0.nodes (list): list of nodes (experiments) in giant component
    """
    # Create a graph using matrix
    G = nx.from_pandas_adjacency(adjacency_matrix)

    # Find giant component
    gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    G0 = G.subgraph(gcc[0])

    # Find nodes of G0 (these will be the index values of the original truth table)
    return list(G0.nodes)

def filter_exp_truth_table(initial_truth_table: pd.DataFrame, 
                           giant_component_nodes_index: list) -> pd.DataFrame:
    """
    Remove rows (experiments) from original truth table that were not in
    the giant component.

    Parameters:
        initial_truth_table (pd.DataFrame): initial truth table with
            rows = experiments and columns = antibodies
        giant_component_nodes_index (list): list of nodes in giant component
    
    Returns:
        filtered_table (pd.DataFrame): truth table with only experiments found
            in the giant component
    """
    # Filter out rows (experiments) from original truth table
    filtered_table = initial_truth_table[initial_truth_table.index.isin(giant_component_nodes_index)]
    return filtered_table

def entire_reference_table(antibodies: list, 
                           idExperiment: list) -> pd.DataFrame:
    """
    Queries the database to retrieve reference CITE-Seq data based on
    the provided antibodies and experiments. 

    Parameters:
        antibodies (list): list of antibodies to query for
        idExperiment (list): list of experiments to query for
    
    Returns:
        ref_df (pd.DataFrame): unformatted table containing normalized
        protein counts for all applicable antibodies and experiment IDs
    """
    warnings.filterwarnings('ignore')
    params = _config()
    conn = mysql.connector.connect(**params)
    cursor = conn.cursor()
    
    ab_placeholders = ','.join(['%s'] * len(antibodies))
    idExperiment_placeholders = ','.join(['%s'] * len(idExperiment))
    
    ab_idExperiment_query = """SELECT cells.idCell, antigen_expression.normValue, 
            antibodies.abTarget, antibodies.idAntibody, cells.idCL, experiments.idExperiment
            FROM cells
            INNER JOIN experiments ON cells.idExperiment = experiments.idExperiment
            INNER JOIN antigen_expression ON cells.idCell = antigen_expression.idCell
            INNER JOIN antibodies ON antigen_expression.idAntibody = antibodies.idAntibody
            WHERE antibodies.idAntibody IN (%s)
                AND cells.idExperiment IN (%s);""" % (ab_placeholders, idExperiment_placeholders)

    parameters = antibodies.copy()
    parameters.extend(idExperiment)

    ref_df = pd.read_sql(sql=ab_idExperiment_query, params=parameters, con=conn)
    
    if ref_df is None or len(ref_df.index) == 0:
        return ("No reference data matching these parameters was found")
         
    return ref_df

def format_reference_table(entire_reference_table: pd.DataFrame, 
                           na_threshold: float = 1) -> pd.DataFrame:
    """
    Takes in unformatted table from database and converts to 
    CITE-Seq format (cells x antibodies), with added "idCL" column
    for the cell type of a particular cell (row). 

    Additional filtering can be applied based on the number of 
    missing values (NA) for a given row (cell).

    Parameters:
        entire_reference_table (pd.DataFrame): unformatted table
            containing all normalized protein values for each antibody
        na_threshold (float): percentage of NAs needed to filter
            out a row from the entire reference table

    Returns:
        filtered_df (pd.DataFrame): formatted reference table where
            rows = cells and columns = antibodies. The final column
            is "idCL", which represents the cell type of the cell.
    """

    antibodies = list(set(entire_reference_table['idAntibody'])) # don't use abTarget here!

    # Pivot the DataFrame with idCellOriginal as index
    pivot_df = entire_reference_table.pivot_table(index=['idCell', 'idCL', 'idExperiment'],
                              columns='idAntibody', # don't use abTarget here!
                              values='normValue',
                              aggfunc='first').reset_index()
    
    # Set idCellOriginal as the index
    pivot_df.set_index('idCell', inplace=True)
    
    # Reorder columns, moving 'idCL' to the end
    df_columns = antibodies.copy()

    # Remove any multilevel index
    pivot_df.columns.name = None
    
    df_columns.append('idCL')
    df_columns.append('idExperiment')
    
    pivot_df = pivot_df[df_columns]
    
    # Set the threshold for NaNs (90%)
    threshold_percentage = na_threshold
    
    # Calculate the percentage of NaNs in specified columns for each row
    nan_percentage = pivot_df[antibodies].isna().sum(axis=1) / len(antibodies)
    
    # Filter rows based on the threshold
    # Filters out cells that have more than %threshold of NaN's in their column values
    filtered_df = pivot_df[nan_percentage < threshold_percentage]
    
    return filtered_df

def downsample(entire_reference_table: pd.DataFrame, 
               size: int = 50) -> pd.DataFrame:
    """
    Downsamples the large reference table to 10,000 rows (cells). 
    Performs downsampling in a controlled randomized order, where
    proportions of each cell type in the original table are 
    kept in the downsampled table. 

    Parameters:
        entire_reference_table (pd.DataFrame): original large table
            containing all cells for the given antibodies
        size (int): the minimum number of cells needed to define 
            a cell type population
    
    Returns:
        entire_reference_table (pd.DataFrame): downsampled 
            reference table
    """

    total_num_cells = len(entire_reference_table.index)
    
    # Downsample if number of rows exceeds 10,000
    if total_num_cells <= 10000:
        return entire_reference_table
    else:
        cells_to_keep = []
        combined_dfs = []
        
        # Find all unique idCLs in the table
        unique_idCLs = list(set(entire_reference_table['idCL']))

        # For each idCL, calculate the number of cells present
        for idCL in unique_idCLs:
            idCL_cells = entire_reference_table.loc[entire_reference_table['idCL'] == idCL]

            # Find number of cells for this cell type
            list_of_cells = list(idCL_cells.index)
            num_idCL_cells = len(list_of_cells)

            # Calculate an adjusted sample_amount for each idCL population to choose from
            sample_amount = ((num_idCL_cells)/(total_num_cells)) * 10000

            # Round up the sample amount
            sample_amount_rounded_up = math.ceil(sample_amount)
            
            if sample_amount_rounded_up > size:
                # Randomly sample this number of cells from this idCL population
                sampled_population_index = sample(list_of_cells, sample_amount_rounded_up)

                # Add these cells to cells_to_keep
                cells_to_keep.extend(sampled_population_index)

                # Create a df for just these cells in this cell type
                temp_df = idCL_cells.loc[sampled_population_index]

                # Add this to combined_dfs
                combined_dfs.append(temp_df)
                
            else:
                # If the sample amount was below our threshold, take 50 of the cells remaining
                if num_idCL_cells > size:
                    # Take 50 of these cells
                    smaller_sampled_population_index = sample(list_of_cells, size)

                    # Add these cells to cells_to_keep
                    cells_to_keep.extend(smaller_sampled_population_index)

                    # Create a df for just these cells in this cell type
                    temp_df = idCL_cells.loc[smaller_sampled_population_index]
                    
                    # Add this to combined_dfs
                    combined_dfs.append(temp_df)
                    
                # If there is not even 50 cells in the population, take whatever remains
                else:
                    remaining_sample_population_index = list_of_cells
                    
                    # Add these cells to cells_to_keep
                    cells_to_keep.extend(remaining_sample_population_index)

                    # Create a df for just these cells in this cell type
                    temp_df = idCL_cells.loc[remaining_sample_population_index]
                    
                    # Add this to combined_dfs
                    combined_dfs.append(temp_df)

        if len(cells_to_keep) > 10000:
            reduced_cells_to_keep = sample(cells_to_keep, 10000)
            return entire_reference_table.loc[pd.Index(reduced_cells_to_keep)]
        else:
            return entire_reference_table.loc[pd.Index(cells_to_keep)]
        
def remove_all_zeros_or_na(protein_df: pd.DataFrame) -> pd.DataFrame:
    """
    Removes any rows (cells) or columns (antibodies) that are all 0s or NAs

    Parameters:
        protein_df (pd.DataFrame): table containing protein count data

    Returns:
        filtered_df (pd.DataFrame): filtered table after removing
            any rows/columns with all 0s or NAs
    """
    # Check if any row in the DataFrame has all NAs or all zeros
    rows_to_exclude = protein_df.apply(lambda row: all(row.isna() | (row == 0)), axis=1)
    
    # Filter out rows to keep only those that do not meet the exclusion conditions
    filtered_df_rows = protein_df[~rows_to_exclude]

    # Check if any column in the DataFrame has all NAs or all zeros
    columns_to_exclude = filtered_df_rows.apply(lambda col: all(col.isna() | (col == 0)), axis=0)
    
    # Filter out columns to keep only those that do not meet the exclusion conditions
    filtered_df = filtered_df_rows.loc[:, ~columns_to_exclude]

    # Display the modified DataFrame
    return filtered_df

def downsample_reference_table(antibody_pairs: list,
                               idBTO: list = None,
                               idExperiment: list = None,
                               parse_option = 1,
                               pairwise_threshold: float = 1.0,
                               na_threshold: float = 1.0,
                               population_size: int = 50) -> pd.DataFrame:
    """
    Wrapper function for generating an appropriate reference CITE-Seq dataset
    used for STvEA annotation transfer

    Parameters:
        antibody_pairs (list): list of pairs, where a pair is:
            [ ["CD3e", "AB_2800722"], 
              ["CD8", "AB_2814271"], 
              ["CD56", "AB_2814309"] ... ]
        idBTO (list): list of tissue IDs to filter for
        idExperiment (list): list of experiment IDs to filter for
        parse_option (int): level of strictness when searching for
            antibodies in the database from an experiment spreadsheet
        pairwise_threshold (float): threshold to decide whether a pairwise distance
            can be used to account for a given experiment
        na_threshold (float): percentage of NAs needed to filter
            out a row from the entire reference table
        population_size (int): the minimum number of cells needed to define 
            a cell type population

    Returns:
        clean_downsampled_df (pd.DataFrame): formatted and downsampled
            reference CITE-Seq dataset where rows = cells, columns = antibodies
    """
    
    # Create antibody dictionary
    antibodies_dict = parse_antibodies(antibody_pairs, 
                                       option = parse_option, 
                                       idBTO=idBTO, 
                                       idExperiment=idExperiment) 

    antibodies = list(antibodies_dict.keys())
                                   
    # Find truth table for all experiments and antibodies
    initial_tb = find_exp_truth_table(antibodies=antibodies, idBTO=idBTO, idExperiment=idExperiment)
               
    # Drop all columns (antibodies) that only contain 0s (false)
    first_ab_drop = drop_antibodies(initial_tb)

    # # Calculate pairwise distances between idExperiment and remaining antibodies
    # pairwise_dist = pairwise_distance_matrix(first_ab_drop)
                                   
    # # Create an adjacency matrix based on filtered pairwise dist values
    # adj_matrix = pairwise_adjacency_matrix(pairwise_dist, pairwise_threshold)
                                   
    # # Create graph and find the giant component's nodes. These nodes values correspond
    # # to the index values found in first_ab_drop
    # gc_nodes_index = find_giant_component(adj_matrix)

    # # Filter out experiments from first_ab_drop that were not in the giant component
    # final_tb = filter_exp_truth_table(first_ab_drop, gc_nodes_index)

    # # Drop all columns again (antibodies) that only contain 0s (false)
    # second_ab_drop = drop_antibodies(final_tb)

    # Use the remaining antibodies and experiments to query for cells in the database
    antibodies_to_query = [ab for ab in first_ab_drop.columns if ab != "idExperiment"] 
    experiments_to_query = list(first_ab_drop["idExperiment"])
                                   
    big_table = entire_reference_table(antibodies_to_query, experiments_to_query)

    # Convert big_table into CITE-SEQ format (cells x antibodies)
    # Last 2 columns: idCL, idExperiment
    format_df = format_reference_table(big_table, na_threshold=na_threshold)
    renamed_format_df = format_df.rename(columns=antibodies_dict, inplace=False)
                                   
    # Downsample to 10k cells
    downsampled_df = downsample(renamed_format_df, size=population_size)

    # Remove any rows or columns that are all 0s or NAs
    clean_downsampled_df = remove_all_zeros_or_na(downsampled_df)
    
    return clean_downsampled_df