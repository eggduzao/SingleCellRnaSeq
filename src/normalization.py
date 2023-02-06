from __future__ import print_function
"""
Normalization Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python

# Internal
from src.util import AuxiliaryFunctions

# External
import scanpy as sc


###################################################################################################
# Goba Class
###################################################################################################

class Normalization():
    """This class represents TODO.

    *Keyword arguments:*

      - argument1 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.

      - argument2 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.
    """

    def __init__(self, anndata_expression_matrix, temporary_folder_name):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Class attributes  
        self.anndata_expression_matrix = anndata_expression_matrix
        self.temporary_folder_name = temporary_folder_name
        self.scanpy_fig_folder_name = "./figures/"


    def normalization_workflow(self, total_count_normalization_target = 1e4, normalization_method = "log1p", variable_genes_min_mean = 0.0125, variable_genes_max_mean = 3,
                               variable_genes_min_disp = 0.5):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """


        # Total Count Normalization
        sc.pp.normalize_total(adata, target_sum = total_count_normalization_target)

        # Normalize by log1p
        if(normalization_method == "log1p")
            sc.pp.log1p(adata)

        # Identify highly variable genes
        sc.pp.highly_variable_genes(adata, min_mean = variable_genes_min_mean, max_mean = variable_genes_max_mean, min_disp = variable_genes_min_disp)

        # Plot of highly variable genes vs others
        sc.pl.highly_variable_genes(adata, show = False, save='_5.pdf')
        AuxiliaryFunctions.save_to_another_folder("violinplot_of_statistics.pdf", self.scanpy_fig_folder_name, self.temporary_folder_name) ###### HERE

        # Filter highly variable genes
        adata = adata[:, adata.var.highly_variable]

        # Regress out effects of total counts of genes and % mito genes
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

        # Scale each gene to unit variance. Clip values exceeding standard deviation 10.
        sc.pp.scale(adata, max_value=10)

































