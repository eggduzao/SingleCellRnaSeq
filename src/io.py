from __future__ import print_function
"""
InputOutput Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import random

# Internal

# External
from scanpy import read_csv, read_h5ad
import pandas as pd


###################################################################################################
# Goba Class
###################################################################################################

class InputOutput():
    """This class represents TODO.

    *Keyword arguments:*

      - argument1 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.

      - argument2 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.
    """

    def __init__(self, input_file_name, input_file_type, temporary_location, temporary_file_type, output_file_type, files_to_remove, gene_alias_instance, seed = None):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        
        # Seed
        if(seed):
            random.seed(seed)
        
        # Class Parameters
        self.input_file_name = input_file_name
        self.input_file_type = input_file_type
        self.temporary_location = temporary_location
        self.temporary_file_type = temporary_file_type
        self.output_file_type = output_file_type
        
        # Auxiliary Parameters
        self.files_to_remove = files_to_remove
        self.gene_alias_instance = gene_alias_instance
        
    ###############################################################################################
    # Reading
    ###############################################################################################

    def read(self, io_input_file_is_reversed = False):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        
        # CSV File
        if(self.input_file_type == "csv"):
        
            # Treat reversed input file
            reversed_input_file_name = os.path.join(self.temporary_location, "reversed_input_file_name.csv")
            if(io_input_file_is_reversed):
                pd.read_csv(self.input_file_name, header=None).T.to_csv(reversed_input_file_name, header=False, index=False)
                self.files_to_remove.append(reversed_input_file_name)
            else:
                reversed_input_file_name = self.input_file_name

            # Treat gene aliases
            alias_input_file_name = os.path.join(self.temporary_location, "alias_input_file_name.csv")
            self.gene_alias_instance.put_gene_names_in_csv_matrix(reversed_input_file_name, alias_input_file_name, io_input_file_is_reversed)
            self.files_to_remove.append(alias_input_file_name)
    
            return read_csv(alias_input_file_name, first_column_names = True)


        # H5AD File
        elif(self.input_file_type == "h5ad"):
        
            return read_h5ad(self.input_file_name)



    ###############################################################################################
    # Writing
    ###############################################################################################





