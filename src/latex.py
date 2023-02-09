from __future__ import print_function
"""
Latex Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os

# Internal

# External


###################################################################################################
# Goba Class
###################################################################################################

class Latex():
    """This class represents TODO.

    *Keyword arguments:*

      - argument1 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.

      - argument2 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.
    """

    def __init__(self, temporary_location, output_report_location):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        
        # Class Parameters
        self.temporary_location = temporary_location
        self.output_report_location = output_report_location
        self.report_file_name = None
        
        
    def open_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        self.report_file_name = open(os.path.join(self.output_report_location, "report.tex"), "w")
    
    def close_report(self): 
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        self.report_file_name.close()


    def add_figure(self, figure_file_name, figure_caption, figure_label, position = "htb", text_width = 1): 
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        self.report_file_name.write("\\begin{figure}["+position+"]\n")
        self.report_file_name.write("\\centering\n")
        self.report_file_name.write("\\includegraphics[width=" + str(text_width) + "\\textwidth]{" + figure_file_name + "}\n")
        self.report_file_name.write("\\centering\n")
        self.report_file_name.write("\\textcolor{black}{\\textbf{\\caption{" + figure_caption + "}}}\\label{" + figure_label + "}\n")
        self.report_file_name.write("\\end{figure}\n\n")

    def create_first_page(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        
        self.report_file_name.write("\\documentclass{article}\n")
        self.report_file_name.write("\\usepackage{reportstylebasic}\n")
        self.report_file_name.write("\\usepackage{subfigure}\n")
        self.report_file_name.write("\\usepackage{longtable}\n")
        self.report_file_name.write("\\usepackage{hyperref} %hyperlinks into use\n\n")
        
        self.report_file_name.write("\\setlength{\\parskip}{1em} %vertical break after paragraph\n")
        self.report_file_name.write("\\setlength\\parindent{0pt} %no \n\n")
        
        self.report_file_name.write("%default font:\n")
        self.report_file_name.write("\\usepackage[scaled=0.92]{helvet}\n")
        self.report_file_name.write("\\renewcommand{\\familydefault}{\\sfdefault}\n\n")

        self.report_file_name.write("%------------------Document----------------------------\n\n")

        self.report_file_name.write("\\begin{document}\n\n")

        self.report_file_name.write("%\\input{newcommands}\n")
        self.report_file_name.write("\\input{./texsettingsandmaterials/texparameters.tex}\n")
        self.report_file_name.write("\\newcommand{\\thedate}{\\today}\n\n")
        
        self.report_file_name.write("%-------------------------TitlePage------------------------\n")
        self.report_file_name.write("%\\begin{titlepage}\n")
        self.report_file_name.write("%\\end{titlepage}\n\n")
        
        self.report_file_name.write("\\thispagestyle{empty} % no page number on cover page\n\n")
        
        self.report_file_name.write("\\section*{Expression profiling study\\\\ project \\project}\n\n")
        
        self.report_file_name.write("\\begin{figure}[htb]\n")
        self.report_file_name.write("    \\centering\n")
        self.report_file_name.write("    \\includegraphics[width=1\\textwidth]{texsettingsandmaterials/figures/cropped-logo_black_wide-2.png}\n")
        self.report_file_name.write("    \\centering\n")
        self.report_file_name.write("    %\\textcolor{Orange}{\\textbf{\\caption{Figure Description}}}\\label{fig:1}\n")        
        self.report_file_name.write("    \\end{figure}\n")
        self.report_file_name.write("\\textbf{\mbcauthors}\n")
        self.report_file_name.write("\\subsection*{Turku Bioscience Centre \\\\ University of Turku and \\AA bo Akademi University \\\\ http://elolab.utu.fi \\\\ mbc@utu.fi}\n\n")        
        
        self.report_file_name.write("\\thedate\n\n")
        
        self.report_file_name.write("\\clearpage\n\n")
        
        self.report_file_name.write("\\tableofcontents\n\n")
        
        self.report_file_name.write("\\clearpage\n\n")


    def create_tools_versions_report(self, tool_version_dictionary):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        
        # Starting session
        self.report_file_name.write("\\section{Tool's Versions}\n\n")

        # Starting table
        self.report_file_name.write("\\begin{table}[htb]\n")
        self.report_file_name.write("\\centering\n")
        self.report_file_name.write("\\textcolor{black}{\\textbf{\\caption{Version of the tools and packages used.}}}\label{tab:tool.versions}\n")
        self.report_file_name.write("\\begin{tabular}{|l|l|l|}\n")
        self.report_file_name.write("\\hline\n")
        
        # Iterating through values
        self.report_file_name.write("\\textbf{Tool} & \\textbf{Version} & \\textbf{Citation} \\\\ \\hline\n")
        for key, value in tool_version_dictionary.items():
            self.report_file_name.write(str(key) + " & " + str(value[0]) + " & " + str(value[1]) + " \\\\ \\hline\n")
            
        # Ending Table
        self.report_file_name.write("\\end{tabular}\n")
        self.report_file_name.write("\\end{table}\n\n")
        
        self.report_file_name.write("\\clearpage\n\n")

    def create_alignment_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        
        self.report_file_name.write("\\section{Samples}\n\n")
        
        self.report_file_name.write("Some Text Here\n\n")
        
        self.report_file_name.write("\\subsubsection*{Result files}\n\n")
        
        self.report_file_name.write("Some more text here\n\n")
        
        self.report_file_name.write("\\clearpage\n\n")
            

    def create_raw_data_preprocessing_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        self.report_file_name.write("\\section{Samples}\n\n")
        
        self.report_file_name.write("Some Text Here\n\n")
        
        self.report_file_name.write("\\subsubsection*{Result files}\n\n")
        
        self.report_file_name.write("Some more text here\n\n")
        
        self.report_file_name.write("\\clearpage\n\n")

    def create_quality_control_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        
        self.report_file_name.write("\\section{Quality Control}\n\n")
        
        # Figure - Top expressed genes before normalization
        figure_file_name = os.path.join(self.temporary_location, "top_expressed_genes.pdf")
        figure_caption = "Placeholder"
        figure_label = "top.expressed.genes"
        self.add_figure(figure_file_name, figure_caption, figure_label, position = "htb", text_width = 0.8)
        
        # Figure - Violin plot of statistics
        figure_file_name = os.path.join(self.temporary_location, "violinplot_of_statistics.pdf")
        figure_caption = "Placeholder"
        figure_label = "violinplot.of.statistics"
        self.add_figure(figure_file_name, figure_caption, figure_label, position = "htb", text_width = 1)
        
        # Figure - Scatterplot of total counts vs mitochondrial counts
        figure_file_name = os.path.join(self.temporary_location, "scatter_total_mito.pdf")
        figure_caption = "Placeholder"
        figure_label = "scatter.total.mito"
        self.add_figure(figure_file_name, figure_caption, figure_label, position = "htb", text_width = 0.5)
        
        # Figure - Scatterplot of total counts vs total genes by counts
        figure_file_name = os.path.join(self.temporary_location, "scatter_total_gbc.pdf")
        figure_caption = "Placeholder"
        figure_label = "scatter.total.gbc"
        self.add_figure(figure_file_name, figure_caption, figure_label, position = "htb", text_width = 0.5)
        
        
        self.report_file_name.write("\\clearpage\n\n")


    def create_normalization_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        
        self.report_file_name.write("\\section{Normalization}\n\n")
        
        # Figure - Top expressed genes before normalization
        figure_file_name = os.path.join(self.temporary_location, "normalized_variable_genes.pdf")
        figure_caption = "Placeholder"
        figure_label = "normalized.variable.genes"
        self.add_figure(figure_file_name, figure_caption, figure_label, position = "htb", text_width = 1)
        
        self.report_file_name.write("\\clearpage\n\n")
        

    def create_feature_selection_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        
        self.report_file_name.write("\\section{Feature Selection}\n\n")
        
        self.report_file_name.write("Placeholder.\n\n")
        
        self.report_file_name.write("\\clearpage\n\n")


    def create_dimensionality_reduction_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        self.report_file_name.write("\\section{Dimensionality Reduction}\n\n")
        
        self.report_file_name.write("Placeholder.\n\n")
        
        self.report_file_name.write("\\clearpage\n\n")
        

    def create_clustering_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        self.report_file_name.write("\\section{Clustering}\n\n")
        
        self.report_file_name.write("Placeholder.\n\n")
        
        self.report_file_name.write("\\clearpage\n\n")
        

    def create_annotation_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        pass

    def create_data_integration_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        pass
        
    def create_de_analysis_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        pass

    def create_compositional_analysis_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        pass

    def create_gsea_analysis_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        pass

    def create_pseudotemporal_ordering_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        pass

    def create_rna_velocity_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        pass

    def create_lineage_tracing_report(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        pass

    def create_bibliography(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        self.report_file_name.write("\\section{File format}\n\n")
        
        self.report_file_name.write("\\begin{itemize}\n")
        self.report_file_name.write(" \\item[•] Point 1 about file format.\n")
        self.report_file_name.write(" \\item[•] Point 2 about file format.\n")
        self.report_file_name.write(" \\item[•] Point 3 about file format.\n")
        self.report_file_name.write("\\end{itemize}\n\n")
        
        self.report_file_name.write("\\section{General notifications}\n")
        self.report_file_name.write("\\begin{itemize}\n")
        self.report_file_name.write(" \\item[•] General Notification 1.\n")
        self.report_file_name.write(" \\item[•] General Notification 2.\n")
        self.report_file_name.write(" \\item[•] General Notification 3.\n")
        self.report_file_name.write("\\end{itemize}\n\n")
        
        self.report_file_name.write("\\clearpage\n\n")
         
        self.report_file_name.write("%Bibliography......................................................................................\n")
        self.report_file_name.write("%\\setlinespacing{1.66}\n")
        self.report_file_name.write("\\bibliographystyle{plain}\n")
        self.report_file_name.write("\\addcontentsline{toc}{section}{\\refname}\n")
        self.report_file_name.write("\\bibliography{texsettingsandmaterials/HTBWorkflow}{}\n")
        self.report_file_name.write("%%.................................................................................................\n\n")
        
        self.report_file_name.write("\\clearpage\n\n")
        
        self.report_file_name.write("\\end{document}\n\n")




