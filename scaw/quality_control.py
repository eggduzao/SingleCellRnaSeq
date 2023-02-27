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
from src.util import AuxiliaryFunctions, InputMatrixColumnType

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

    def __init__(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Class attributes
        self.scanpy_fig_folder_name = "./figures/"
        

    def quality_control_workflow(self, anndata_expression_matrix, temporary_folder_name, top_expressed_genes = 20, min_genes_per_cell = 200,
                                 min_cells_per_gene = 3, max_number_of_genes_by_counts = 2500, max_percentage_mithocondrial_counts = 5):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        
        # Make var names unique in case `var_names='gene_ids'` in `sc.read_10x_mtx`
        anndata_expression_matrix.var_names_make_unique()

        # Plot highly expressed genes
        sc.pl.highest_expr_genes(anndata_expression_matrix, n_top = top_expressed_genes, show = False, save='_top_expressed_genes.pdf')
        AuxiliaryFunctions.save_to_another_folder("top_expressed_genes.pdf", self.scanpy_fig_folder_name, temporary_folder_name)

        # Basic filtering
        sc.pp.filter_cells(anndata_expression_matrix, min_genes = min_genes_per_cell)
        sc.pp.filter_genes(anndata_expression_matrix, min_cells = min_cells_per_gene)

        # Annotate mithocontrial genes and perform mithocontrial QC
        anndata_expression_matrix.var['mt'] = anndata_expression_matrix.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(anndata_expression_matrix, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        # Violin plot of statistics
        sc.pl.violin(anndata_expression_matrix, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show = False, save='_violinplot_of_statistics.pdf')
        AuxiliaryFunctions.save_to_another_folder("violinplot_of_statistics.pdf", self.scanpy_fig_folder_name, temporary_folder_name)

        # Scatter plot of total counts vs percentage of mithocondrial counts
        sc.pl.scatter(anndata_expression_matrix, x='total_counts', y='pct_counts_mt', show = False, save='_scatter_total_mito.pdf')
        AuxiliaryFunctions.save_to_another_folder("scatter_total_mito.pdf", self.scanpy_fig_folder_name, temporary_folder_name)
        
        # Scatter plot of total counts vs number of genes by counts
        sc.pl.scatter(anndata_expression_matrix, x='total_counts', y='n_genes_by_counts', show = False, save='_scatter_total_gbc.pdf')
        AuxiliaryFunctions.save_to_another_folder("scatter_total_gbc.pdf", self.scanpy_fig_folder_name, temporary_folder_name)

        # Filter by slicing the AnnData object
        anndata_expression_matrix = anndata_expression_matrix[anndata_expression_matrix.obs.n_genes_by_counts < max_number_of_genes_by_counts, :]
        anndata_expression_matrix = anndata_expression_matrix[anndata_expression_matrix.obs.pct_counts_mt < max_percentage_mithocondrial_counts, :]


