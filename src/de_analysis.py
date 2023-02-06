from __future__ import print_function
"""
DE Analysis Module
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
import pandas as pd

###################################################################################################
# Goba Class
###################################################################################################

class DEAnalysis():
    """This class represents TODO.

    *Keyword arguments:*

      - argument1 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.

      - argument2 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.
    """

    def __init__(self, placeholder):
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


    def finding_marker_genes_workflow(self, clustering_method = "leiden", statistical_method = "wilcox", number_of_genes_to_plot = 25,
                                      top_genes_for_table = 10):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Finding marker genes based on a t-test
        if(statistical_method == "ttest"):
            sc.tl.rank_genes_groups(adata, clustering_method, method = 't-test')
            sc.pl.rank_genes_groups(adata, n_genes = number_of_genes_to_plot, sharey = False, show = False, save = "_top_marker_genes.pdf")
            AuxiliaryFunctions.save_to_another_folder("top_marker_genes.pdf", self.scanpy_fig_folder_name, self.temporary_folder_name)

        # Finding marker genes based on a Wilcoxon-Mann-Whitney test
        elif(statistical_method == "wilcox"):
            sc.tl.rank_genes_groups(adata, clustering_method, method = 'wilcoxon')
            sc.pl.rank_genes_groups(adata, n_genes = number_of_genes_to_plot, sharey = False, show = False, save = "_top_marker_genes.pdf")
            AuxiliaryFunctions.save_to_another_folder("top_marker_genes.pdf", self.scanpy_fig_folder_name, self.temporary_folder_name)

        # Finding marker genes based on a logistic regression
        elif(statistical_method == "logreg"):
            sc.tl.rank_genes_groups(adata, clustering_method, method = 'logreg')
            sc.pl.rank_genes_groups(adata, n_genes = number_of_genes_to_plot, sharey = False, show = False, save = "_top_marker_genes.pdf")
            AuxiliaryFunctions.save_to_another_folder("top_marker_genes.pdf", self.scanpy_fig_folder_name, self.temporary_folder_name)


        # Show the top X ranked genes in a table
        top_marker_genes_table = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(top_genes_for_table)
        top_marker_genes_table.to_csv(os.path.join(self.temporary_folder_name, "top_marker_genes_table.csv"))

        """
# Show the top X ranked genes in a table with scores
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(5)

# Violin plot for cluster vs rest
sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)

# Violin plot of list of genes by cluster
sc.pl.violin(adata, ['CST3', 'NKG7', 'PPBP'], groupby='leiden')

# Compare to each of the single clusters
#sc.tl.rank_genes_groups(adata, 'leiden', groups=['0'], reference='1', method='wilcoxon')
#sc.pl.rank_genes_groups(adata, groups=['0'], n_genes=20)

# Violin plot for each of the single clusters
#sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)

        """
