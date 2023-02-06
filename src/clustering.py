from __future__ import print_function
"""
Clustering Module
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

class Clustering():
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

    def clustering_workflow(self, clustering_method = "leiden", list_of_genes_of_interest = []):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Clustering method
        if(clustering_method == "leiden"):
            sc.tl.leiden(self.anndata_expression_matrix)
        else:
            sc.tl.louvain(self.anndata_expression_matrix)

        # Plot the clusters
        sc.pl.umap(self.anndata_expression_matrix, color=[clustering_method] + list_of_genes_of_interest, show = False, save="_clustering_umap_data.pdf")
        AuxiliaryFunctions.save_to_another_folder("clustering_umap_data.pdf", self.scanpy_fig_folder_name, self.temporary_folder_name)
























