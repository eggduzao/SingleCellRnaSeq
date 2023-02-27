from __future__ import print_function
"""
Arguments Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
import optparse
import subprocess
import multiprocessing

# Internal
from scaw.__version__ import __version__
from scaw.util import PassThroughOptionParser, AuxiliaryFunctions, InputMatrixColumnType

# External

###################################################################################################
# ArgumentParser Class
###################################################################################################

class ArgumentParser():
    """This class represents TODO.

    *Keyword arguments:*

      - argument1 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.

      - argument2 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.
    """

    def __init__(self, error_handler):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """
    
        # Initializing class objects
        self.error_handler = error_handler
        self.usage_message = None
        self.version_message = None
        self.parser = None
        self.options = None
        self.arguments = None

        # Load options and arguments
        self.load_usage_message()
        self.load_version_message()
        self.load_parser()
        self.load_options()
        self.load_options_and_arguments()
        self.option_argument_validity_check()

    #############################################################################
    # Main Operations
    #############################################################################

    def load_usage_message(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Usage message
        self.usage_message = ("%prog [options] <Input Matrix File> <Temporary Location> <Output Location>\n\n"

                              "-----------------------------------------------------------------------------\n"
                              "Scaw is a workflow designed to analyze single-cell RNA-seq data\n"
                              "from their source and generate a report with statistics.\n\n"

                              "The <Input Matrix File> can be in one of these formats:\n"
                              "- Format 1.\n" # TODO
                              "  Format 2.\n" # TODO
                              "- Format 3.\n\n" # TODO

                              "Scaw's documentation can be found at:\n"
                              "WEBSITE\n\n" # TODO

                              "For more information, please refer to:\n"
                              "WEBSITE\n\n" # TODO

                              "For further questions or comments please refer to:\n"
                              "ORIGINAL PAPER\n"  # TODO
                              "-----------------------------------------------------------------------------")


    def load_version_message(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Version message
        self.version_message = "SCAW - Single-Cell Analysis Workflow. Version: " + str(__version__)

    def load_parser(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Initializing Option Parser
        self.parser = PassThroughOptionParser(usage = self.usage_message, version = self.version_message)

    def add_option(self, option_alias, option_variable, option_type, meta_type, default_option, help_message):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Creating option alias
        option_alias = "-" + option_alias

        # Creating option name based on variable name
        option_name = "--" + option_variable.replace("_", "-")

        # Correct meta_type to uppercase
        meta_type = meta_type.upper()

        # Hidden options from help message
        if(not help_message): help_message = optparse.SUPPRESS_HELP
    
        # Add option
        if(option_type == "bool" and default_option):
            self.parser.add_option(option_alias, option_name, dest = option_variable, action="store_false",
                                   default=default_option, help=help_message)
        elif(option_type == "bool"):
            self.parser.add_option(option_alias, option_name, dest = option_variable, action="store_true",
                                   default=default_option, help=help_message)
        else:
            self.parser.add_option(option_alias, option_name, dest = option_variable, type = option_type,
                                   metavar = meta_type, default=default_option, help=help_message)

    def load_options(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
        
        *Return:*
    
          - return -- A return.
        """

        # Options

        organism_help = ("The organism in which the analysis is being applied to. It can be one of "
                         "the following: hg19, hg38, mm9, mm10, zv9 or zv10. If you are analyzing other "
                         "organism,please check the manual for the necessary files and configuration.")
        self.add_option("o", "organism", "string", "STRING", "hg19", organism_help)

        input_type_help = ("The format of the input file. SCAW currently supports: XXXXXXXXXXXXXXXX ")
        self.add_option("i", "input_type", "string", "STRING", None, input_type_help)

        cores_help = ("Number of cores the user would like to dedicate for the parallel steps of "
                      "SCAW pipeline. Default to half the total number of CPU cores.")
        self.add_option("c", "cores", "int", "INT", multiprocessing.cpu_count() / 2, cores_help)

        output_type_help = ("The desired format of the output file. Since only one matrix is generated by "
                            "SCAW, possible alternatives are: XXXXXXXXXXXXXXXXXXXXXX.")
        self.add_option("y", "output_type", "string", "STRING", "sparse", output_type_help)

        seed_help = ("Select the same pseudo-random number generator seed to be able to"
                     "replicate an experiment (same output is guaranteed for the same seed).")
        self.add_option("s", "seed", "int", "INT", 123, seed_help)

        # Hidden Options

        io_input_file_is_reversed_help = None
        self.add_option("A", "io_input_file_is_reversed", "bool", "BOOLFALSE", False, io_input_file_is_reversed_help)

        pre_remove_threshold_help = None
        self.add_option("B", "pre_remove_threshold", "float", "FLOAT", 0.0, pre_remove_threshold_help)

        sica_pvalue_threshold_help = None
        self.add_option("C", "sica_pvalue_threshold", "float", "FLOAT", 0.99, sica_pvalue_threshold_help)

        sica_bottom_bin_ext_range_help = None
        self.add_option("D", "sica_bottom_bin_ext_range", "string", "STRING", "3,10", sica_bottom_bin_ext_range_help)

        """
        # Examples:
    
        parameter_1_help = ("Parameter 1 does a lot of things like bla bla bla bla bla bla bla bla bl"
                            "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
        self.add_option("a", "parameter1", "string", "PATH", os.getcwd(), parameter_1_help)

        parameter_2_help = ("Parameter 2 does a lot of things like bla bla bla bla bla bla bla bla bl"
                            "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
        self.add_option("b", "parameter_2", "string", "string", None, parameter_2_help)

        parameter_3_help = ("Parameter 3 does a lot of things like bla bla bla bla bla bla bla bla bl"
                            "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
        self.add_option("c", "parameter3", "int", "INT", 1, parameter_3_help)
    
        parameter_4_help = ("Parameter 4 does a lot of things like bla bla bla bla bla bla bla bla bl"
                        "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
        self.add_option("d", "parameter4_44_4", "float", "FLOAT", 0.1, parameter_4_help)

        parameter_5_help = ("Parameter 5 does a lot of things like bla bla bla bla bla bla bla bla bl"
                            "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
        self.add_option("e", "parameter5", "bool", "BOOLFALSE", None, parameter_5_help)

        parameter_6_help = None
        self.add_option("parameter6", "string", "STRING", "STRING", parameter_6_help)
    """

    def load_options_and_arguments(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Processing Options
        self.options, self.arguments = self.parser.parse_args()

    def option_argument_validity_check(self):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Verify number of arguments
        if(len(sys.argv) == 1):
            self.parser.print_help(sys.stderr)
            sys.exit(1)
        if(len(self.arguments) != 3):
            pass
            # self.error_handler.throw_error("TODO") # TODO

        # Arguments
        input_matrix = self.arguments[0]
        output_matrix = self.arguments[1]
        output_matrix_directory = os.path.dirname(os.path.abspath(os.path.expanduser(output_matrix)))
        output_contacts = self.arguments[2]
        output_contacts_directory = os.path.dirname(os.path.abspath(os.path.expanduser(output_contacts)))

        # Verify arguments
        # TODO

        # Options
        organism = self.options.organism
        input_type = self.options.input_type
        cores = self.options.cores
        output_type = self.options.output_type
        seed = self.options.seed

        # Verify operational options
        # TODO


    def read_input_manifest_file(self):
        """Reads the input manifest format table. The table has the following format:

    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Appending input table
        input_table = []
        input_file_name = self.arguments[0]
        input_file = open(input_file_name, "r")
        for line in input_file:
        
            # Comment and header lines
            if(len(line) < 5 or line[0] == "#"): continue
            if("GroupID" in line): continue
            
            # Canonical Lines
            ll = line.strip().split("\t")
            ll[InputMatrixColumnType.INPUTFILE] = os.path.abspath(os.path.expanduser(ll[InputMatrixColumnType.INPUTFILE]))
            
            # Temporary condition folder
            temporary_condition_folder = os.path.abspath(os.path.expanduser(os.path.join(self.arguments[1], ll[InputMatrixColumnType.CONDITION])))
            if(not os.path.exists(temporary_condition_folder)):
                os.makedirs(temporary_condition_folder)
            ll.append(temporary_condition_folder)
            
            # Output condition folder
            output_condition_folder = os.path.abspath(os.path.expanduser(os.path.join(self.arguments[2], ll[InputMatrixColumnType.CONDITION])))
            if(not os.path.exists(temporary_condition_folder)):
                os.makedirs(temporary_condition_folder)
            ll.append(output_condition_folder)
            
            input_table.append(ll)
        input_file.close()

        # Returning objects
        return input_table


    #############################################################################
    # Auxiliary Operations
    #############################################################################

    def create_temporary_directory(self, input_file_name, temporary_location):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Creating temorary directory as: temporary folder / input file name
        input_contact_matrix_name = os.path.splitext(os.path.basename(input_file_name))[0]
        temporary_directory = None
        try:
            temporary_directory = os.path.join(temporary_location, input_contact_matrix_name)
            temp_creation_command = ["mkdir", "-p", temporary_directory]
            temp_creation_process = subprocess.run(temp_creation_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
        except Exception:
            raise
            # self.error_handler.throw_error("TODO") # TODO

        # Returning the name of the temporary directory
        return temporary_directory

    def create_output_directory(self, output_location):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Creating output directory
        try:
            output_creation_command = ["mkdir", "-p", output_location]
            output_creation_process = subprocess.run(output_creation_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
        except Exception:
            raise
            # self.error_handler.throw_error("TODO") # TODO


