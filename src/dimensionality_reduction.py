from __future__ import print_function
"""
Dimensionality reduction Module
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

class DimensionalityReduction():
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


    def dimensionality_reduction_workflow(self, number_of_neighbors = 10, number_of_pcs = 40, list_of_genes_of_interest = [],
                                          clustering_method = "leiden"):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        
        # Computing the neighborhood graph
        sc.pp.neighbors(self.anndata_expression_matrix, n_neighbors = number_of_neighbors, n_pcs = number_of_pcs)

        # Embedding the neighborhood graph
        if(clustering_method == "leiden"):
            sc.tl.leiden(self.anndata_expression_matrix)
        else:
            sc.tl.louvain(self.anndata_expression_matrix)
        sc.tl.paga(self.anndata_expression_matrix)
        sc.pl.paga(self.anndata_expression_matrix, plot = False)  # remove `plot=False` if you want to see the coarse-grained graph
        sc.tl.umap(self.anndata_expression_matrix, init_pos = "paga")

        # Compute UMAP
        sc.tl.umap(self.anndata_expression_matrix)

        # Plot UMAP with raw data
        sc.pl.umap(self.anndata_expression_matrix, color = list_of_genes_of_interest, show = False, save="_UMAP_raw_data.pdf")
        AuxiliaryFunctions.save_to_another_folder("UMAP_raw_data.pdf", self.scanpy_fig_folder_name, self.temporary_folder_name)

        # Plot UMAP with normalized data
        sc.pl.umap(self.anndata_expression_matrix, color = list_of_genes_of_interest, use_raw=False, show = False, save="_UMAP_normalized_data.pdf")
        AuxiliaryFunctions.save_to_another_folder("UMAP_normalized_data.pdf", self.scanpy_fig_folder_name, self.temporary_folder_name)


