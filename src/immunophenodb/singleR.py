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
        deltas_df (pd.DataFrame): dataframe containing certainty values for 
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
            srLabels = sr.SingleR(test=r_df, ref=HPCA_rds, labels=dollar(HPCA_rds, "label.ont"))

    # Else, keep current gene names
    elif num_rows != num_ensembl:
        print("No ensembl IDs found...")
        
        with (robjects.default_converter + pandas2ri.converter).context():
            print("Converting pandas dataframe into R dataframe...")
            r_df = robjects.conversion.get_conversion().py2rpy(rna_counts)
            print("Running SingleR...")
            srLabels = sr.SingleR(test=r_df, ref=HPCA_rds, labels=dollar(HPCA_rds, "label.ont"))

    # Convert SingleR R data frame to a pandas df
    with (robjects.default_converter + pandas2ri.converter).context():
        print("Converting R dataframe back into pandas dataframe...")
        pd_df = robjects.conversion.get_conversion().rpy2py(r('as.data.frame')(srLabels))
        pd_labels_delta_df = pd_df[['labels', 'delta.next']]

        labels_df = pd.DataFrame(pd_labels_delta_df.loc[:, 'labels'])
        deltas_df = pd.DataFrame(pd_labels_delta_df.loc[:, 'delta.next'])

    print("Finished running SingleR.")
    return labels_df, deltas_df

def annotate_cells(IPD,
                   ref_ver: int = 75):
    """
    Annotates cells using RNA data from an ImmunoPhenoData object
    Parameters:
        IPD (ImmunoPhenoData object): object that must contain gene data along
            with protein data
        ref_ver (int): version number for reference dataset
    """
    # If RNA data is provided, we can run singleR to get the idCLs and certainty metrics
    if IPD.gene_cleaned is not None:
        rna_counts = IPD.gene_cleaned

        # Setup singleR
        setup_singleR()

        # Run singleR
        labels_df, deltas_df = singleR(rna_counts, ref_ver)

        IPD.raw_cell_labels = labels_df
        IPD.label_certainties = deltas_df
    else:
        raise Exception("Error. No RNA data found. Unable to annotate cells.")
    
    return