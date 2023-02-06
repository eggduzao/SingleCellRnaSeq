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
from src.quality_control import QualityControl
from src.normalization import Normalization
from src.feature_selection import FeatureSelection
from src.dimensionality_reduction import DimensionalityReduction
from src.clustering import Clustering
from src.de_analysis import DEAnalysis

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
import seaborn as se
import matplotlib as mpl
import matplotlib.pyplot as plt

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
    
    # Scanpy initializations
    sc.settings.verbosity = 0
    sc.settings.set_figure_params(dpi=80, facecolor='white')
    
    ###############################################################################################
    # Parameters
    ###############################################################################################
    
    # Initializing argument parser
    argument_parser = ArgumentParser(error_handler)
    args = argument_parser.arguments
    opts = argument_parser.options

    # Arguments
    input_matrix_file_name = os.path.abspath(os.path.expanduser(args[0]))
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
    
    print(io_input_file_is_reversed)
    
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
    
    temporary_file_type = "h5ad"
    
    # Reading input file
    inputoutput_instance = InputOutput(input_matrix_file_name, input_type, temporary_location, temporary_file_type, output_type, files_to_remove, gene_alias, seed)
    
    # Reading file
    anndata_expression_matrix = None
    anndata_expression_matrix = inputoutput_instance.read(io_input_file_is_reversed = io_input_file_is_reversed)
    
    print(anndata_expression_matrix)
    print(anndata_expression_matrix.var[:10])
    
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
    seaborn_version = se.__version__
    matplotlib_version = mpl.__version__
    
    # Further tool's versions
    # TODO
    
    ###############################################################################################
    # Alignment
    ###############################################################################################    

    # TODO
    
    # Alignment preprocessing time
    alignment_preprocessing_timestamp = time.time()  
    
    
    ###############################################################################################
    # Raw Data Preprocessing
    ###############################################################################################
    
    # TODO
    
    # Raw data preprocessing time
    raw_data_preprocessing_timestamp = time.time()
    
    ###############################################################################################
    # Quality Control
    ###############################################################################################
    
    top_expressed_genes = 20
    min_genes_per_cell = 200
    min_cells_per_gene = 3
    max_number_of_genes_by_counts = 2500
    max_percentage_mithocondrial_counts = 5
    
    # Perform quality control
    quality_control_instance = QualityControl(anndata_expression_matrix, temporary_location)
    quality_control_instance.quality_control_workflow(top_expressed_genes = top_expressed_genes,
                                                      min_genes_per_cell = min_genes_per_cell,
                                                      min_cells_per_gene = min_cells_per_gene,
                                                      max_number_of_genes_by_counts = max_number_of_genes_by_counts,
                                                      max_percentage_mithocondrial_counts = max_percentage_mithocondrial_counts)
    
    
    # Quality control time
    quality_control_timestamp = time.time()

    
    ###############################################################################################
    # Normalization
    ###############################################################################################
    
    total_count_normalization_target = 1e4
    normalization_method = "log1p"
    variable_genes_min_mean = 0.0125
    variable_genes_max_mean = 3
    variable_genes_min_disp = 0.5
    clip_values_max_std = 10


    # Perform normalization
    normalization_instance = Normalization(anndata_expression_matrix, temporary_location)
    normalization_instance.normalization_workflow(total_count_normalization_target = total_count_normalization_target,
                                                  normalization_method = normalization_method,
                                                  variable_genes_min_mean = variable_genes_min_mean,
                                                  variable_genes_max_mean = variable_genes_max_mean,
                                                  variable_genes_min_disp = variable_genes_min_disp,
                                                  clip_values_max_std = clip_values_max_std)
    
    # Normalization
    normalization_timestamp = time.time()
    
    
    ###############################################################################################
    # Feature Selection
    ###############################################################################################
    
    list_of_genes_of_interest = ["CST3", "MAFP5"]
    log_variance_ratio = True

    # Perform feature selection
    feature_selection_instance = FeatureSelection(anndata_expression_matrix, temporary_location)
    feature_selection_instance.feature_selection_workflow(list_of_genes_of_interest = list_of_genes_of_interest, 
                                                          log_variance_ratio = log_variance_ratio)

    # Feature selection time
    feature_selection_timestamp = time.time()
    
    
    ###############################################################################################
    # Dimensionality Reduction
    ###############################################################################################
    
    number_of_neighbors = 10
    number_of_pcs = 40
    list_of_genes_of_interest = ["CST3", "NKG7"]
    clustering_method = "leiden"

    # Perform dimensionality reduction
    dimensionality_reduction_instance = DimensionalityReduction(anndata_expression_matrix, temporary_location)
    dimensionality_reduction_instance.dimensionality_reduction_workflow(number_of_neighbors = number_of_neighbors,
                                                                        number_of_pcs = number_of_pcs,
                                                                        list_of_genes_of_interest = list_of_genes_of_interest,
                                                                        clustering_method = clustering_method)

    # Dimensionality reduction time
    dimensionality_reduction_timestamp = time.time()
    
    
    ###############################################################################################
    # Clustering
    ###############################################################################################

    clustering_method = "leiden"
    list_of_genes_of_interest = ["CST3", "NKG7"]    

    # Perform clustering
    clustering_instance = Clustering(anndata_expression_matrix, temporary_location)
    clustering_instance.clustering_workflow(clustering_method = clustering_method, 
                                            list_of_genes_of_interest = list_of_genes_of_interest)

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
    
    # Performing DE Analysis
    de_analysis_instance = DEAnalysis(anndata_expression_matrix, temporary_location)

    # Finding marker genes
    de_analysis_instance.finding_marker_genes_workflow(clustering_method = "leiden",
                                                       statistical_method = "wilcox",
                                                       number_of_genes_to_plot = 25,
                                                       top_genes_for_table = 10)
    
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
    
    
   
