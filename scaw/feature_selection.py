from __future__ import print_function
"""
Feature selection Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python

# Internal
from scaw.util import AuxiliaryFunctions

# External
import scanpy as sc


###################################################################################################
# Goba Class
###################################################################################################

class FeatureSelection():
    """This class represents TODO.

    *Keyword arguments:*

      - argument1 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.

      - argument2 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.
    """

    def __init__(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Class attributes
        self.scanpy_fig_folder_name = "./figures/"


    def feature_selection_workflow(self, anndata_expression_matrix, temporary_folder_name, list_of_genes_of_interest = [], log_variance_ratio = True):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Select features with PCA
        sc.tl.pca(anndata_expression_matrix, svd_solver = "arpack")

        # Make scatterplot with list of interesting genes
        for gene in list_of_genes_of_interest:
            try:
                sc.pl.pca(anndata_expression_matrix, color = gene, show = False, save = "_PCA_expression_"+gene+".pdf")
                AuxiliaryFunctions.save_to_another_folder("PCA_expression_"+gene+".pdf", self.scanpy_fig_folder_name, temporary_folder_name)
            except Exception:
              # Error - gene not found
              pass

        # Plot the variance ratio of the scatterplot
        sc.pl.pca_variance_ratio(anndata_expression_matrix, log = log_variance_ratio, show = False, save="_PCA_variance_ratio.pdf")
        AuxiliaryFunctions.save_to_another_folder("PCA_variance_ratio.pdf", self.scanpy_fig_folder_name, temporary_folder_name)











