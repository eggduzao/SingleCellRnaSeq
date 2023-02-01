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
import random

# Internal

# External
from scanpy import read_csv, read_h5ad



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

    def __init__(self, input_file_name, input_file_type, temporary_file_type, output_file_type, seed = None):
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
        self.temporary_file_type = temporary_file_type
        self.output_file_type = output_file_type
        
        
        """
        # Operation type
        if(operation_type == "read"):
        
            # Reading csv file
            if(input_file_type = "csv"):
                self.read_csv(input_file_name)
                
            # Reading h5ad file
            elif(input_file_type == "h5ad"):
                self.read_h5ad(input_file_name)
        
        elif(operation_type == "write"):
            pass
        """

        
    ###############################################################################################
    # Reading
    ###############################################################################################

    def read_csv(self, input_file_name, delimiter = ",", first_column_names = None, dtype = "float32", transpose = False):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        if(transpose):
            # ssss
            pass
        else:
            read_csv(input_file_name, delimiter, first_column_names, dtype)


    def read_h5ad(self, input_file_name, backed = None, as_sparse = ()):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
        read_h5ad(input_file_name, backed, as_sparse)



    ###############################################################################################
    # Writing
    ###############################################################################################





