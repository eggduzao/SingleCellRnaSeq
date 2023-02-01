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
from src.util import ErrorHandler, JuicerCommand, CoolerCommand, ChromosomeSizes
from src.arguments import ArgumentParser

# External


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
    pre_min_contig_removed_bins = opts.pre_min_contig_removed_bins # TOREMOVE
    pre_remove_threshold = opts.pre_remove_threshold # TOREMOVE
    sica_pvalue_threshold = opts.sica_pvalue_threshold # TOREMOVE
    sica_bottom_bin_ext_range = [int(e) for e in opts.sica_bottom_bin_ext_range.split(",")] # TOREMOVE


    # Setting seed
    random.seed(seed)
    
    # Config data classes
    juicer_command = JuicerCommand() # TOREMOVE
    cooler_command = CoolerCommand() # TOREMOVE
    chromosome_sizes = ChromosomeSizes(organism)
    
    ###############################################################################################
    # Input
    ###############################################################################################
    
    
    # TODO
    
    
    ###############################################################################################
    # Tool's Versions
    ###############################################################################################
    
    
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
    
    
   
