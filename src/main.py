from __future__ import print_function
"""
Main Module
===================
The main module which will call all the functions of the workflow in order.

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
from collections import OrderedDict

# Internal
from src.__version__ import __version__
from src.util import ErrorHandler, JuicerCommand, CoolerCommand, GeneAlias, ReportConfiguration
from src.arguments import ArgumentParser
from src.io import InputOutput
from src.quality_control import QualityControl
from src.normalization import Normalization
from src.feature_selection import FeatureSelection
from src.dimensionality_reduction import DimensionalityReduction
from src.clustering import Clustering
from src.de_analysis import DEAnalysis
from src.latex import Latex

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
    """Entry point to the execution of SCAW.
    
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
    input_matrix_file_name = argument_parser.read_input_manifest_file()
    general_temporary_location = os.path.abspath(os.path.expanduser(args[1]))
    general_output_location = os.path.abspath(os.path.expanduser(args[2]))

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
    report_configuration = ReportConfiguration()
    
    ###############################################################################################
    # Input
    ###############################################################################################
    
    temporary_file_type = "h5ad"
    
    # Reading input file
    inputoutput_instance = InputOutput(input_matrix_file_name, input_type, general_temporary_location, temporary_file_type, output_type, files_to_remove, gene_alias, seed)
    
    # Reading file
    anndata_expression_matrix_vector = None
    anndata_expression_matrix_vector = inputoutput_instance.read(io_input_file_is_reversed = io_input_file_is_reversed)
    
    # Sample information
    sample_information_matrix = []
    counter = 0
    for anndata_expression_matrix in anndata_expression_matrix_vector:
        condition_name = input_matrix_file_name[counter][2] # Condition name
        number_of_cells = len(anndata_expression_matrix.obs) # Number of cells
        number_of_genes = len(anndata_expression_matrix.var) # Number of genes
        sample_information_matrix.append([condition_name, number_of_cells, number_of_genes]) # [Number of cells, number of genes]
        counter += 1

    
    ###############################################################################################
    # Tool's Versions
    ###############################################################################################
    
    # Python version
    python_version = sys.version.split(" ")[0]
    
    # Tool's version
    tool_version = __version__
    
    # Python packages versions
    numpy_version = np.__version__
    scipy_version = si.__version__
    pandas_version = pd.__version__
    scanpy_version = sc.__version__
    anndata_version = an.__version__
    umap_version = um.__version__
    sklearn_version = sk.__version__
    statsmodels_version = st.__version__
    igraph_version = ig.__version__
    louvain_version = lo.__version__
    pynndescent_version = nn.__version__
    seaborn_version = se.__version__
    matplotlib_version = mpl.__version__
    
    # Further tool's versions
    # TODO
    
    # Dictionary of versions
    tool_version_dictionary = OrderedDict()
    tool_version_dictionary["Python"] = [python_version, "~\\cite{van1995python}"]
    tool_version_dictionary["SCAW"] = [tool_version, "~"]
    tool_version_dictionary["Numpy"] = [numpy_version, "~\\cite{harris2020numpy}"]
    tool_version_dictionary["Scipy"] = [scipy_version, "~\\cite{virtanen2020scipy}"]
    tool_version_dictionary["Pandas"] = [pandas_version, "~\\cite{mckinney2010pandas}"]
    tool_version_dictionary["Scanpy"] = [scanpy_version, "~\\cite{wolf2018scanpy}"]
    tool_version_dictionary["AnnData"] = [anndata_version, "~\\cite{virshup2021anndata}"]
    tool_version_dictionary["UMAP"] = [umap_version, "~\\cite{mcinnes2018umap}"]
    tool_version_dictionary["Scikit-Learn"] = [sklearn_version, "~\\cite{pedregosa2011sklearn}"]
    tool_version_dictionary["Statsmodel"] = [statsmodels_version, "~\\cite{seabold2010statsmodels}"]
    tool_version_dictionary["IGraph"] = [igraph_version, "~\\cite{csardi2006igraph}"]
    tool_version_dictionary["Louvain"] = [louvain_version, "~\\cite{aynaud2020louvain}"]
    tool_version_dictionary["PyNNDescent"] = [pynndescent_version, "~\\cite{wei2006PyNNDescent}"]
    tool_version_dictionary["Seaborn"] = [seaborn_version, "~\\cite{waskom2021seaborn}"]
    tool_version_dictionary["Matplotlib"] = [matplotlib_version, "~\\cite{hunter2007matplotlib}"]

    
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

    """

    qc_top_expressed_genes = 20
    qc_min_genes_per_cell = 200
    qc_min_cells_per_gene = 3
    qc_max_number_of_genes_by_counts = 2500
    qc_max_percentage_mithocondrial_counts = 5
    
    # Perform quality control
    quality_control_instance = QualityControl(anndata_expression_matrix, temporary_location)
    quality_control_instance.quality_control_workflow(top_expressed_genes = qc_top_expressed_genes,
                                                      min_genes_per_cell = qc_min_genes_per_cell,
                                                      min_cells_per_gene = qc_min_cells_per_gene,
                                                      max_number_of_genes_by_counts = qc_max_number_of_genes_by_counts,
                                                      max_percentage_mithocondrial_counts = qc_max_percentage_mithocondrial_counts)
    
    
    # Quality control time
    quality_control_timestamp = time.time()

    
    ###############################################################################################
    # Normalization
    ###############################################################################################
    
    nm_total_count_normalization_target = 1e4
    nm_normalization_method = "log1p"
    nm_variable_genes_min_mean = 0.0125
    nm_variable_genes_max_mean = 3
    nm_variable_genes_min_disp = 0.5
    nm_clip_values_max_std = 10


    # Perform normalization
    normalization_instance = Normalization(anndata_expression_matrix, temporary_location)
    normalization_instance.normalization_workflow(total_count_normalization_target = nm_total_count_normalization_target,
                                                  normalization_method = nm_normalization_method,
                                                  variable_genes_min_mean = nm_variable_genes_min_mean,
                                                  variable_genes_max_mean = nm_variable_genes_max_mean,
                                                  variable_genes_min_disp = nm_variable_genes_min_disp,
                                                  clip_values_max_std = nm_clip_values_max_std)
    
    # Normalization
    normalization_timestamp = time.time()
    

    ###############################################################################################
    # Feature Selection
    ###############################################################################################
    
    #fs_list_of_genes_of_interest = ["MFAP5", "MAFP5"]
    fs_list_of_genes_of_interest = ["SAT1", "FTL"]
    fs_log_variance_ratio = True

    # Perform feature selection
    feature_selection_instance = FeatureSelection(anndata_expression_matrix, temporary_location)
    feature_selection_instance.feature_selection_workflow(list_of_genes_of_interest = fs_list_of_genes_of_interest, 
                                                          log_variance_ratio = fs_log_variance_ratio)

    # Feature selection time
    feature_selection_timestamp = time.time()
    
    
    ###############################################################################################
    # Dimensionality Reduction
    ###############################################################################################
    
    dr_number_of_neighbors = 10
    dr_number_of_pcs = 40
    #dr_list_of_genes_of_interest = ["MFAP5", "SPARCL1"]
    dr_list_of_genes_of_interest = ["SAT1", "FTL"]
    dr_clustering_method = "leiden"

    # Perform dimensionality reduction
    dimensionality_reduction_instance = DimensionalityReduction(anndata_expression_matrix, temporary_location)
    dimensionality_reduction_instance.dimensionality_reduction_workflow(number_of_neighbors = dr_number_of_neighbors,
                                                                        number_of_pcs = dr_number_of_pcs,
                                                                        list_of_genes_of_interest = dr_list_of_genes_of_interest,
                                                                        clustering_method = dr_clustering_method)

    # Dimensionality reduction time
    dimensionality_reduction_timestamp = time.time()
    

    ###############################################################################################
    # Clustering
    ###############################################################################################

    cl_clustering_method = "leiden"
    #cl_list_of_genes_of_interest = ["MFAP5", "SPARCL1"] 
    cl_list_of_genes_of_interest = ["SAT1", "FTL"]    

    # Perform clustering
    clustering_instance = Clustering(anndata_expression_matrix, temporary_location)
    clustering_instance.clustering_workflow(clustering_method = cl_clustering_method, 
                                            list_of_genes_of_interest = cl_list_of_genes_of_interest)

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
    
    clustering_method = "leiden"
    statistical_method = "wilcox"
    number_of_genes_to_plot = 25
    top_genes_for_table = 10
    groups_of_interest = ["0", "5"]
    #list_of_genes_of_interest = ["MFAP5", "SPARCL1"]
    list_of_genes_of_interest = ["SAT1", "FTL"]
    
    # Performing DE Analysis
    de_analysis_instance = DEAnalysis(anndata_expression_matrix, temporary_location)

    # Finding marker genes
    de_analysis_instance.finding_marker_genes_workflow(clustering_method = clustering_method,
                                                       statistical_method = statistical_method,
                                                       number_of_genes_to_plot = number_of_genes_to_plot,
                                                       top_genes_for_table = top_genes_for_table,
                                                       groups_of_interest = groups_of_interest,
                                                       list_of_genes_of_interest = list_of_genes_of_interest)
    
    # DE Analysis time
    de_analysis_timestamp = time.time()
    
    """ 
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
    # Create Report
    ###############################################################################################

    # Report variables
    output_report_location = os.path.join(general_output_location, "report")
    report_configuration.move_latex_files(output_report_location)
    
    # LaTeX report
    latex_instance = Latex(general_temporary_location, output_report_location)
    latex_instance.open_report()
    latex_instance.create_first_page()
    latex_instance.create_tools_versions_report(tool_version_dictionary)
    latex_instance.create_sample_information_report(sample_information_matrix)
    #latex_instance.create_alignment_report()
    #latex_instance.create_raw_data_preprocessing_report()
    #latex_instance.create_quality_control_report()
    #latex_instance.create_normalization_report()
    #latex_instance.create_feature_selection_report()
    #latex_instance.create_dimensionality_reduction_report()
    #latex_instance.create_clustering_report()
    #latex_instance.create_annotation_report()
    #latex_instance.create_data_integration_report()
    #latex_instance.create_de_analysis_report()
    #latex_instance.create_compositional_analysis_report()
    #latex_instance.create_gsea_analysis_report()
    #latex_instance.create_pseudotemporal_ordering_report()
    #latex_instance.create_rna_velocity_report()
    #latex_instance.create_lineage_tracing_report()
    latex_instance.create_bibliography()
    latex_instance.close_report()
    
    # Running LaTeX report
    report_configuration.run_report(output_report_location)

    
    ###############################################################################################
    # Deleting objects
    ###############################################################################################

    
    # Files to remove
    for filename in files_to_remove:
        command = "rm -rf " + filename
        os.system(command)
    

    print("Completed in :" + str(lineage_tracing_timestamp - alignment_preprocessing_timestamp) + "s")

    
 
    
    
   
