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

SCI_CRUNCH_BASE = "http://www.scicrunch.org"
SCI_RRID_ENDPOINT = "/resolver/"
SCI_FILE = ".json"

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
    idCL_dict['not_idCLs'] = []
    
    # Loop over each antibody_id (key) in the outer dictionary
    for antibody_id, results in lmm_results.items():
        # results is: LMM table, aggregate_df, proportion of background cells not in idCLs list
        lmm_table = results[0]
        aggregate_df = results[1]
        bg_proportion = results[2]
        
        cellType_coeff.append(float(lmm_table.loc['cellType']['Coef.']))
        stdError.append(float(lmm_table.loc['cellType']['Std.Err.']))
        p_values.append(float(lmm_table.loc['cellType']['P>|z|']))

        for key in parents:
            value = aggregate_df.loc[key]['proportion']
            idCL_dict[key].append(value)

        idCL_dict['not_idCLs'].append(bg_proportion)
            
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
    final_df.sort_values('coeff', inplace=True)

    return final_df

def get_antibodies(node_fam_dict: dict, 
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
    # Collapse dict into one list of cell types
    parents = list(node_fam_dict.keys())
    children = [item for sublist in list(node_fam_dict.values()) for item in sublist]
    family = list(set(parents + children))
    
    # Retrieve unique antibodies from idCLs
    antibodies = get_antibodies_idCL(family, idBTO, idExperiment)

    # Run Spark job to on antibodies
    ab_lmm_output = run_spark_query_ab_lmm(antibodies, node_fam_dict)

    if len(ab_lmm_output) == 0:
        return ("No antibodies found. Please retry with fewer cell types or different parameters.")
    
    # Convert output into a readable dataframe
    result_df = lmm_ab_list_to_df(ab_lmm_output, parents)

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
        idCL_dict["background"] = []

        for idCL, table_and_proportion_tuple in idCL_lmm_result.items():
            df = table_and_proportion_tuple[0]
            proportion = table_and_proportion_tuple[1]
            
            cellType_coeff.append(float(df.loc['cellType']['Coef.']))
            stdError.append(float(df.loc['cellType']['Std.Err.']))
            p_values.append(float(df.loc['cellType']['P>|z|']))

            idCL_dict["background"].append(proportion)
            
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
    
        # Sort by smallest coefficient first
        final_df.sort_values('coeff', inplace=True)

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