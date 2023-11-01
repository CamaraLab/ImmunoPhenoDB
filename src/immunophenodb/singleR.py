import pandas as pd
import logging
from importlib.resources import files

import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr, data
from rpy2.robjects import pandas2ri, r
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger

import gzip
import io
import requests
from pyensembl import EnsemblRelease
from gtfparse import read_gtf

import requests
import multiprocessing
import multiprocess

def setup_singleR():
    """
    Installs necessary packages for running SingleR via rpy2
    """
    # Choosing a CRAN Mirror
    utils = importr('utils')
    utils.chooseCRANmirror(ind=70)

    # Install required packages if not already 
    package_names = ("BiocManager", "SummarizedExperiment", "SingleR")
    if all(rpackages.isinstalled(x) for x in package_names):
        print("All R packages installed...")
        have_package = True
    else:
        print("Missing R packages...")
        have_package = False

    if not have_package:
        rpy2_logger.setLevel(logging.ERROR)
        print("Installing missing R packages...")
        robjects.r('install.packages("BiocManager", quietly=TRUE)')
        robjects.r('BiocManager::install("SummarizedExperiment", ask=FALSE, quietly=TRUE)')
        robjects.r('BiocManager::install("SingleR", ask=FALSE, quietly=TRUE)')

def get_ensdb_ref(ref_ver: int) -> pd.DataFrame:
    """
    Retrieves a reference dataset from the Human Primary Cell Atlas database

    Parameters:
        ref_ver (int): version number of the reference dataset
    
    Returns:
        ref_genes_final (pd.DataFrame): dataframe containing reference data
    """
    # Retrieve the download link containing the reference data
    ref_link = EnsemblRelease(ref_ver).gtf_url

    print(f"Downloading reference dataset for version:{ref_ver}...")
    # Save gz data into memory
    response = requests.get(ref_link, stream=True)  
    response_gz_file = response.content
    f = io.BytesIO(response_gz_file)
    print("Download Complete...")

    # Extract using gzip into another memory container, then open as pandas df
    with gzip.GzipFile(fileobj=f) as fh:
        ref_df = read_gtf(fh)
    
    # Extract the gene_name (SYMBOL) and gene_id columns, removing duplicates
    ref_genes = ref_df[['gene_name', 'gene_id']].drop_duplicates()

    # Set second column (gene ids) to be the index
    ref_genes_final = ref_genes.set_index(ref_genes.columns[1])

    return ref_genes_final

def singleR(rna_count_df: pd.DataFrame, 
            ref_ver: int = None) -> pd.DataFrame:
    """
    Runs SingleR using rpy2

    Parameters:
        rna_count_df (pd.DataFrame): dataframe containing RNA data
        ref_ver (int): version number for reference dataset
    
    Returns:
        labels_df (pd.DataFrame): dataframe containing cell annotations as
            cell ontologies (EMBL-EBI conventions, CL:XXXXXXX)
        singleR_df (pd.DataFrame): dataframe containing certainty values for 
            each cell annotation
    """
    # Instead of installing Human Primary Cell Atlas from celldex, load in RDS file (SummarizedExperiment object)
    HPCA_rds_path = str(files('immunophenodb.data').joinpath('HPCA.rds'))
    readRDS = robjects.r['readRDS']
    HPCA_rds = readRDS(file=HPCA_rds_path)

    # Use $ operator from R
    base = importr('base')
    dollar = base.__dict__["$"]

    # Import SingleR 
    sr = importr('SingleR')

    # Read RNA as pandas df
    print("Loading RNA counts into a pandas dataframe...")
    rna_counts = rna_count_df.T

    # CHECK: Do ensembl gene IDs exist in the RNA count?
    print("Checking for ensembl IDs...")
    num_rows = len(rna_counts.index)
    num_ensembl = rna_counts.index.str.startswith('ENSG').sum()

    # Convert all ensembl gene IDs to their common names
    if num_rows == num_ensembl:
        print("Ensembl IDs found. Converting into common names...")
        # Get reference gene IDs 
        ref_genes = get_ensdb_ref(ref_ver)

        # Get corresponding common gene names
        common_name_index = ref_genes.loc[rna_counts.index]['gene_name']

        # Replace old index with common names
        rna_counts_sub = rna_counts.set_index(common_name_index)
        with (robjects.default_converter + pandas2ri.converter).context():
            print("Converting pandas dataframe into R dataframe...")
            # Convert new pandas df into R dataframe
            r_df = robjects.conversion.get_conversion().py2rpy(rna_counts_sub)
            print("Running SingleR...")
            # Run SingleR
            srLabels = sr.SingleR(test=r_df, ref=HPCA_rds, labels=dollar(HPCA_rds, "label.ont"), prune=False)

    # Else, keep current gene names
    elif num_rows != num_ensembl:
        print("No ensembl IDs found...")
        # Convert pandas df into R data frame
        
        with (robjects.default_converter + pandas2ri.converter).context():
            print("Converting pandas dataframe into R dataframe...")
            # Convert new pandas df into R dataframe
            r_df = robjects.conversion.get_conversion().py2rpy(rna_counts)
            print("Running SingleR...")
            # Run SingleR
            srLabels = sr.SingleR(test=r_df, ref=HPCA_rds, labels=dollar(HPCA_rds, "label.ont"), prune=False)

    # Convert SingleR R data frame to a pandas df
    with (robjects.default_converter + pandas2ri.converter).context():
        print("Converting R dataframe back into pandas dataframe...")
        pd_df = robjects.conversion.get_conversion().rpy2py(r('as.data.frame')(srLabels))

        labels_df = pd.DataFrame(pd_df.loc[:, 'labels'])
        singleR_df = pd.DataFrame(pd_df)

    return labels_df, singleR_df

def convert_idCL_readable(idCL:str) -> str:
    """
    Converts a cell ontology id (CL:XXXXXXX) into a readable cell type name

    Parameters:
        idCL (str): cell ontology ID

    Returns:
        cellType (str): readable cell type name
        
    """
    idCL_params = {
        'q': idCL,
        'exact': 'true',
        'ontology': 'cl',
        'fieldList': 'label',
        'rows': 1,
        'start': 0
    }
    
    res = requests.get("https://www.ebi.ac.uk/ols/api/search", params=idCL_params)
    res_JSON = res.json()
    cellType = res_JSON['response']['docs'][0]['label']
    
    return cellType

def ebi_idCL_map(labels_df: pd.DataFrame) -> dict:
    """
    Converts a list of cell ontology IDs into readable cell type names
    as a dictionary

    Parameters:
        labels_df (pd.DataFrame): dataframe with cell labels from singleR
    
    Returns:
        idCL_map (dict) : dictionary mapping cell ontology ID to cell type
    
    """
    idCL_map = {}
    
    idCLs = set(labels_df["labels"])
    
    for idCL in idCLs:
        idCL_map[idCL] = convert_idCL_readable(idCL)
    
    return idCL_map

def generate_range(size: int, 
                   interval_size: int = 2000,
                   start: int = 0) -> list:
    """
    Splits a range of numbers into partitions based on a
    provided interval size. Partitions will only contain
    the begin and end of a range. Ex: [[0, 10], [10, 20]]
    
    Parameters:
        size (int): total size of the range
        interval_size (int): desired partition sizes 
        start (int): index to begin partitioning the range
    
    Returns:
        start_stop_ranges (list): list of lists, containing
            the start and end of a partition. 
    """
    
    if interval_size <= 1:
        raise ValueError("Chunk or partition size must be greater than 1.")
    
    iterable_range = [i for i in range(size + 1)]

    end = max(iterable_range)
    step = interval_size
    
    start_stop_ranges = []
    
    for i in range(start, end, step): 
        x = i 
        start_stop_ranges.append([iterable_range[x:x+step][0], 
                                  iterable_range[x:x+step+1][-1]])
    
    return start_stop_ranges

def clean_chunks(chunk_list: list) -> list:
    """
    Merges any chunks in a list that only contain one element with 
    the previous chunk. This is because SingleR requires at least 
    two cells at once to run.
    
    Parameters:
        chunk_list (list): list of lists containing start and end values
        
    Returns:
        chunk_list (list): updated list of lists without any lists
            that contain a single element
    """
    final_chunk = chunk_list[-1]
    size_final_chunk = final_chunk[1] - final_chunk[0]
    
    if size_final_chunk == 1:
        # Update the second to last chunk with the last chunk
        # "Merge" by adding the last chunk to the second to last chunk
        # and delete the last chunk
        chunk_list[-2][1] = chunk_list[-2][1] + 1
        chunk_list.pop()
        
    return chunk_list

def singleR_parallel(rna, chunk_size=20000, partition_size=2000):
    """
    Executes SingleR using parallelized batch processing with 
    a provided RNA dataframe.
    
    Large RNA dataframes will be divided into large 'chunks', which 
    are divided into smaller 'partitions'.
    
    Chunks will be evaluated sequentially, while partitions inside a chunk
    will be evaluated in parallel using SingleR.
    
    Parameters:
        rna (pd.DataFrame): sparse pandas dataframe containing RNA data
            with rows (cells) x columns (genes)
        chunk_size (int): size of a chunk within the RNA dataframe
        partition_size (int): size of a partition within a chunk
        
    Returns:
        all_labels (pd.DataFrame): singleR output with cell labels
        all_singleR (pd.DataFrame): singleR output with label certainties
    """
    chunk_ranges = clean_chunks(generate_range(len(rna.index), 
                                               start=0, 
                                               interval_size=chunk_size))
    all_labels_list = []
    all_singleR_list = []
    
    def process_section(section):
        # Check if all columns in section are sparse. If so, convert to dense
        all_types = rna.iloc[section[0] : section[1]].dtypes

        if all(isinstance(t, pd.core.arrays.sparse.dtype.SparseDtype) for t in all_types):
            cells_batch = rna.iloc[section[0] : section[1]].sparse.to_dense()
        # Otherwise, use default dataframe which is not sparse
        else:
            cells_batch = rna.iloc[section[0] : section[1]]

        print(f"Running SingleR on partition section: {section}\n")
        labels_df, singleR_df = singleR(cells_batch)
        return labels_df, singleR_df
    
    for chunk in chunk_ranges:
        print("Chunk section:", chunk)
        partition_ranges = clean_chunks(generate_range(chunk[1], 
                                                       interval_size=partition_size, 
                                                       start=chunk[0]))

        # Create a multiprocessing pool
        with multiprocess.Pool(processes=multiprocessing.cpu_count()) as pool:
            results = pool.map(process_section, partition_ranges)

        for labels_df, singleR_df in results:
            all_labels_list.append(labels_df)
            all_singleR_list.append(singleR_df)
            
    return pd.concat(all_labels_list), pd.concat(all_singleR_list)

def annotate_cells(IPD,
                   chunk_size: int = 20000,
                   partition_size: int = 2000,
                   ref_ver: int = 75):
    """
    Annotates cells using RNA data from an ImmunoPhenoData object
    
    Parameters:
        IPD (ImmunoPhenoData object): object that must contain gene data along
            with protein data
        chunk_size (int): size of a chunk used to divide a large RNA dataframe
        partition_size (int): size of a partition within a chunk used
            when running SingleR in parallel
        ref_ver (int): version number for reference dataset
    """
    
    if IPD._singleR_rna is not None:
        rna_counts = IPD._singleR_rna
        # Setup singleR
        setup_singleR()
        
        # Run parallelized version of SingleR
        labels_df, singleR_df = singleR_parallel(rna_counts,chunk_size, partition_size)
        
        # Create mapping dictionary of idCLs to cell type names
        complete_map_dict = ebi_idCL_map(labels_df)
        
        labels_df['celltype'] = labels_df['labels'].map(complete_map_dict)
        
        # Store in IPD object
        IPD.raw_cell_labels = labels_df
        IPD._temp_labels = labels_df
        
        IPD.label_certainties = singleR_df
        IPD._temp_certainties = singleR_df
        
        print("Finished running SingleR.")
    else:
        raise Exception("Error. No RNA data found. Unable to annotate cells.")
    
    return