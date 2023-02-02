from __future__ import print_function
"""
Quality Control Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python

# Internal
from src.Util import AuxiliaryFunctions

# External
import scanpy as sc


###################################################################################################
# Goba Class
###################################################################################################

class QualityControl():
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
        
        self.anndata_expression_matrix = anndata_expression_matrix
        self.temporary_folder_name = temporary_folder_name
        self.scanpy_fig_folder_name = "./figures/"
        

    def quality_control_workflow(self, top_expressed_genes = 20, min_genes_per_cell = 200, min_cells_per_gene = 3, max_number_of_genes_by_counts = 2500,
                                 max_percentage_mithocondrial_counts = 5):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        
        # Make var names unique in case `var_names='gene_ids'` in `sc.read_10x_mtx`
        self.anndata_expression_matrix.var_names_make_unique()

        # Plot highly expressed genes
        sc.pl.highest_expr_genes(self.anndata_expression_matrix, n_top = top_expressed_genes, show = False, save='_top_expressed_genes.pdf')
        #AuxiliaryFunctions.save_to_another_folder("top_expressed_genes.pdf", self.scanpy_fig_folder_name, self.temporary_folder_name)

        # Basic filtering
        sc.pp.filter_cells(self.anndata_expression_matrix, min_genes = min_genes_per_cell)
        sc.pp.filter_genes(self.anndata_expression_matrix, min_cells = min_cells_per_gene)

        # Annotate mithocontrial genes and perform mithocontrial QC
        self.anndata_expression_matrix.var['mt'] = self.anndata_expression_matrix.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(self.anndata_expression_matrix, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        # Violin plot of statistics
        sc.pl.violin(self.anndata_expression_matrix, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show = False, save='_violinplot_of_statistics.pdf')
        # TODO - Save plot to temporary folder

        # Scatter plot of total counts vs percentage of mithocondrial counts
        sc.pl.scatter(self.anndata_expression_matrix, x='total_counts', y='pct_counts_mt', show = False, save='_scatter_total_mito.pdf')
        # TODO - Save plot to temporary folder
        
        # Scatter plot of total counts vs number of genes by counts
        sc.pl.scatter(self.anndata_expression_matrix, x='total_counts', y='n_genes_by_counts', show = False, save='_scatter_total_gbc.pdf')
        # TODO - Save plot to temporary folder

        # Filter by slicing the AnnData object
        self.anndata_expression_matrix = self.anndata_expression_matrix[self.anndata_expression_matrix.obs.n_genes_by_counts < max_number_of_genes_by_counts, :]
        self.anndata_expression_matrix = self.anndata_expression_matrix[self.anndata_expression_matrix.obs.pct_counts_mt < max_percentage_mithocondrial_counts, :]

    































