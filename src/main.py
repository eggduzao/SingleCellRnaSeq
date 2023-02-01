from __future__ import print_function
"""
Main Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
import time
import warnings
import random

# Internal
from src.__version__ import __version__
from src.util import ErrorHandler, JuicerCommand, CoolerCommand, GeneAlias
from src.arguments import ArgumentParser
from src.io import InputOutput

# External
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as an
import umap as um
import scipy as si
import sklearn as sk
import statsmodels as st
import igraph as ig
import louvain as lo
import pynndescent as nn

###################################################################################################
# Tool Execution
###################################################################################################

def main():
    """Entry point to the execution of Bloom.
    
    *Keyword arguments:*
    
      - None.
    
    *Return:*
    
      - None.
    """

    ###############################################################################
    # Error and Warning Handler + Initializations
    ###############################################################################

    # Initializing ErrorHandler
    error_handler = ErrorHandler()

    # Do not show system warnings
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
    
    # Workflow starting time
    start_timestamp = time.time()
    
    # Files to delete at the end of the workflow
    files_to_remove = []
    
    ###############################################################################################
    # Parameters
    ###############################################################################################
    
    # Initializing argument parser
    argument_parser = ArgumentParser(error_handler)
    args = argument_parser.arguments
    opts = argument_parser.options

    # Arguments
    input_matrix = os.path.abspath(os.path.expanduser(args[0]))
    temporary_location = os.path.abspath(os.path.expanduser(args[1]))
    output_location = os.path.abspath(os.path.expanduser(args[2]))

    # Options
    organism = opts.organism
    input_type = opts.input_type
    cores = opts.cores
    output_type = opts.output_type
    seed = opts.seed

    # Hidden Options
    io_input_file_is_reversed = opts.io_input_file_is_reversed
    
    #pre_min_contig_removed_bins = opts.pre_min_contig_removed_bins # TOREMOVE
    #pre_remove_threshold = opts.pre_remove_threshold # TOREMOVE
    #sica_pvalue_threshold = opts.sica_pvalue_threshold # TOREMOVE
    #sica_bottom_bin_ext_range = [int(e) for e in opts.sica_bottom_bin_ext_range.split(",")] # TOREMOVE


    # Setting seed
    random.seed(seed)
    
    # Config data classes
    gene_alias = GeneAlias(organism)
    juicer_command = JuicerCommand() # TOREMOVE
    cooler_command = CoolerCommand() # TOREMOVE
    
    ###############################################################################################
    # Input
    ###############################################################################################
    
    
    # Reading input file
    temporary_file_type = "h5ad"
    inputoutput_instance = InputOutput(input_matrix, input_type, temporary_file_type, output_type, seed)
    
    # Treat reversed input file
    reversed_input_file_name = os.path.join(temporary_location, "reversed_input_file_name.csv")
    if(io_input_file_is_reversed):
        pd.read_csv(input_file_name, header=None).T.to_csv(reversed_input_file_name, header=False, index=False)
        files_to_remove.append(reversed_input_file_name)
    else:
        reversed_input_file_name = input_matrix

    # Treat gene aliases
    alias_input_file_name = os.path.join(temporary_location, "alias_input_file_name.csv")
    gene_alias.put_gene_names_in_csv_matrix(reversed_input_file_name, alias_input_file_name)
    files_to_remove.append(alias_input_file_name)
    
    ###############################################################################################
    # Tool's Versions
    ###############################################################################################
    
    # Python version
    python_version = sys.version.split(" ")[0]
    
    # Tool's version
    tool_version = __version__
    
    # Python packages versions
    numpy_version = np.__version__
    pandas_version = pd.__version__
    scanpy_version = sc.__version__
    anndata_version = an.__version__
    umap_version = um.__version__
    scipy_version = si.__version__
    sklearn_version = sk.__version__
    statsmodels_version = st.__version__
    igraph_version = ig.__version__
    louvain_version = lo.__version__
    pynndescent_version = nn.__version__
    
    print(python_version)
    print(tool_version)
    
    print(numpy_version)
    print(pandas_version)
    print(scanpy_version)
    print(anndata_version)
    print(umap_version)
    print(scipy_version)
    print(sklearn_version)
    print(statsmodels_version)
    print(igraph_version)
    print(louvain_version)
    print(pynndescent_version)
    
    
    # Further tool's versions
    # TODO
    
    ###############################################################################################
    # Raw Data Preprocessing
    ###############################################################################################
    
    # TODO
    
    # Raw data preprocessing time
    raw_data_preprocessing_timestamp = time.time()
    
    ###############################################################################################
    # Quality Control
    ###############################################################################################
    

    
    # Quality control time
    quality_control_timestamp = time.time()

    
    ###############################################################################################
    # Normalization
    ###############################################################################################
    
    
    # Normalization
    normalization_timestamp = time.time()
    
    
    ###############################################################################################
    # Feature Selection
    ###############################################################################################
    
    # TODO
    
    # Feature selection time
    feature_selection_timestamp = time.time()
    
    
    ###############################################################################################
    # Dimensionality Reduction
    ###############################################################################################
    
    # TODO
    
    # Dimensionality reduction time
    dimensionality_reduction_timestamp = time.time()
    
    
    ###############################################################################################
    # Clustering
    ###############################################################################################
    
    # TODO
    
    # Clustering time
    clustering_timestamp = time.time()
    
    
    ###############################################################################################
    # Annotation
    ###############################################################################################
    
    # TODO
    
    # Annotation time
    annotation_timestamp = time.time()
    
    
    ###############################################################################################
    # Data Integration
    ###############################################################################################
    
    # TODO
    
    # Data integration time
    data_integration_timestamp = time.time()
    
    
    ###############################################################################################
    # DE Analysis
    ###############################################################################################
    
    # TODO
    
    # DE Analysis time
    de_analysis_timestamp = time.time()
    
    
    ###############################################################################################
    # Compositional Analysis
    ###############################################################################################
    
    # TODO
    
    # Compositional analysis
    compositional_analysis_timestamp = time.time()
    
    
    ###############################################################################################
    # GSEA Analysis
    ###############################################################################################
    
    # TODO
    
    # GSEA analysis time
    gsea_analysis_timestamp = time.time()
    
    
    ###############################################################################################
    # Pseudotemporal Ordering
    ###############################################################################################
    
    # TODO
    
    # Pseudotemporal ordering time
    Pseudotemporal_ordering_timestamp = time.time()
    
    
    ###############################################################################################
    # RNA Velocity
    ###############################################################################################
    
    # TODO
    
    # RNA Velocity
    rna_velocity_timestamp = time.time()
    
    
    ###############################################################################################
    # Lineage Tracing
    ###############################################################################################
    
    # TODO
    
    # Lineage tracing time
    lineage_tracing_timestamp = time.time()
    
    
    ###############################################################################################
    # Placeholder
    ###############################################################################################
    
    
   
