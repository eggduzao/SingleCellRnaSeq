from __future__ import print_function
"""
Util Module
===================
The Util classes contains many utilities needed by other classes such as the paths to input files.

Authors: Eduardo G. Gusmao. Partially based on the Util class from RGT (https://www.regulatory-genomics.org/).

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
import glob
import codecs
import optparse
import traceback
import subprocess
import configparser

# Internal

# External
import numpy

###################################################################################################
# Configuration File Handling
###################################################################################################

class ConfigurationFile:
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

        # Fetching scaw_data_path
        self.scaw_data_path = os.path.expanduser(os.getenv("scaw_data", default = os.path.join(os.getenv("HOME"), "scaw_data")))

        # Reading config file directory
        self.scaw_config_file_name = os.path.join(self.scaw_data_path, "data.config")

        # Parsing config file
        self.config = configparser.ConfigParser()
        self.config.read_file(codecs.open(self.scaw_config_file_name, "rU", "utf8"))

        # Reading data directory
        #self.data_dir = os.path.split(data_config_file_name)[0]

class GenomeData(ConfigurationFile):
    """Represent genomic data. Inherits ConfigurationFile."""

    def __init__(self, organism):
        """Initializes GenomeData.

        *Keyword arguments:*

            - organism -- Organism alias.
        """
        ConfigurationFile.__init__(self)
        self.organism = organism
        self.genome = os.path.join(self.data_dir, self.config.get(organism, 'genome'))
        self.chromosome_sizes = os.path.join(self.data_dir, self.config.get(organism, 'chromosome_sizes'))
        self.gene_regions = os.path.join(self.data_dir, self.config.get(organism, 'gene_regions'))
        self.annotation = os.path.join(self.data_dir, self.config.get(organism, 'annotation'))
        self.annotation_dump_dir = os.path.dirname(os.path.join(self.data_dir, self.annotation))
        self.gene_alias = os.path.join(self.data_dir, self.config.get(organism, 'gene_alias'))
        if organism in ["hg19", "hg38", "mm9"]:
            self.repeat_maskers = os.path.join(self.data_dir, self.config.get(organism, 'repeat_maskers'))
        else:
            self.repeat_maskers = None

    def get_organism(self):
        """Returns the current organism."""
        return self.organism

    def get_genome(self):
        """Returns the current path to the genome fasta file."""
        return self.genome

    def get_chromosome_sizes(self):
        """Returns the current path to the chromosome sizes text file."""
        return self.chromosome_sizes

    def get_gene_regions(self):
        """Returns the current path to the gene_regions BED file."""
        return self.gene_regions

    def get_annotation(self):
        """
        Returns the current path to the gencode annotation gtf file.
        """
        return self.annotation

    def get_annotation_dump_dir(self):
        """Returns the current path to the gencode annotation gtf file."""
        return self.annotation_dump_dir

    def get_gene_alias(self):
        """Returns the current path to the gene alias txt file."""
        return self.gene_alias

    def get_repeat_maskers(self):
        """Returns the current path to directory for repeat maskers."""
        if self.repeat_maskers:
            return self.repeat_maskers
        else:
            print("*** There is no repeat masker data for " + self.organism)


class JuicerCommand(ConfigurationFile):
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

        # Configuration file initialization
        ConfigurationFile.__init__(self)
        self.juicer_command = self.config.get("Juicer", "command")
        self.juicer_options = self.config.get("Juicer", "options")
        self.juicer_jar_location = os.path.join(self.scaw_data_path, self.config.get("Juicer", "jar"))

class CoolerCommand(ConfigurationFile):
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

        # Configuration file initialization
        ConfigurationFile.__init__(self)
        self.cooler_command = self.config.get("Cooler", "command")

class ChromosomeSizes(ConfigurationFile):
    """This class represents TODO.

    *Keyword arguments:*

      - argument1 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.

      - argument2 -- Short description. This argument represents a long description. It can be:
        - Possibility 1: A possibility 1.
        - Possibility 2: A possibility 2.
    """

    def __init__(self, organism):
        """Returns TODO.
    
        *Keyword arguments:*
    
          - argument -- An argument.
    
        *Return:*
    
          - return -- A return.
        """

        # Configuration file initialization
        ConfigurationFile.__init__(self)
        self.organism = organism
        self.chromosome_sizes_file_name = os.path.join(self.scaw_data_path, self.config.get("ChromosomeSizes", organism))

        # Creating chromosome sizes dictionary and chromosome list
        self.chromosome_sizes_dictionary = dict()
        chrom_sizes_file = codecs.open(self.chromosome_sizes_file_name, "rU", "utf8")
        for line in chrom_sizes_file:
            ll = line.strip().split("\t")
            self.chromosome_sizes_dictionary[ll[0]] = int(ll[1])
        chrom_sizes_file.close()
        self.chromosome_sizes_list = sorted(self.chromosome_sizes_dictionary.keys())


###################################################################################################
# Argument Parsing
###################################################################################################

class HelpfulOptionParser(optparse.OptionParser):
    """This class represents an OptionParser that prints full help on errors. Inherits OptionParser.

    *Keyword arguments:*

      - OptionParser -- Inherited OptionParser object.
    """

    def error(self, msg):
        """Error handling.
    
        *Keyword arguments:*
    
          - msg -- String containing the error message.
    
        *Return:*
    
          - return -- An error message to the user.
        """
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


class PassThroughOptionParser(HelpfulOptionParser):
    """When unknown arguments are encountered, bundle with largs and try again, until rargs is depleted.
       sys.exit(status) will still be called if a known argument is passed incorrectly (e.g. missing arguments or
       bad argument types, etc.). Inherits HelpfulOptionParser.

    *Keyword arguments:*

      - HelpfulOptionParser -- Inherited HelpfulOptionParser object.
    """

    def _process_args(self, largs, rargs, values):
        """Overwrites _process_args to achieve desired effect.
    
        *Keyword arguments:*
    
          - largs -- The tool's optional arguments.
          - rargs -- The tool's required arguments.
        """
        while rargs:
            try:
                HelpfulOptionParser._process_args(self, largs, rargs, values)
            except (optparse.BadOptionError, optparse.AmbiguousOptionError):
                pass
                # largs.append(err.opt_str)


###################################################################################################
# Warnings and Errors Handling
###################################################################################################

class ErrorHandler:
    """Handles errors in a standardized way.

        *Error Dictionary Standard:*

            Each entry consists of a key+list in the form X:[Y,Z,W] where:

            - X -- The key representing the internal error name.
            - Y -- Error number.
            - Z -- Exit status.
            - W -- Error message to be print.

        *Warning Dictionary Standard:*

            Each entry consists of a key+list in the form X:[Y,Z] where:

            - X -- The key representing the internal warning name.
            - Y -- Warning number.
            - Z -- Warning message to be print.
    """

    def __init__(self):

        self.program_name = os.path.basename(sys.argv[0])

        self.error_dictionary = {
            "DEFAULT_ERROR": [0, 0, "Undefined error. Program terminated with exit status 0."],
            "MOTIF_ANALYSIS_OPTION_ERROR": [1, 0,
                                            "You must define one specific analysis. Run '" + self.program_name + " -h' for help."],

            "FP_WRONG_ARGUMENT": [2, 0,
                                  "You must provide at least one and no more than one experimental matrix as input argument."],
            "FP_WRONG_EXPMAT": [3, 0,
                                "The experimental matrix could not be loaded. Check if it is correctly formatted and that your python version is >= 2.7."],
            "FP_ONE_REGION": [4, 0, "You must provide one 'regions' bed file in the experiment matrix."],
            "FP_NO_DNASE": [5, 0, "You must provide one 'reads' file termed DNASE or ATAC in the experiment matrix."],
            "FP_NO_HISTONE": [6, 0,
                              "You must provide at least one 'reads' file not termed DNASE or ATAC (histone modification) in the experiment matrix."],
            "FP_NB_HMMS": [7, 0,
                           "You must provide one HMM file or X HMM files where X is the number of histone tracks detected in the experiment matrix."],
            "FP_HMM_FILES": [8, 0,
                             "Your HMM file could not be read. If you did not provide any HMM file, you may not have installed scikit correctly."],
            "FP_BB_CREATION": [9, 0,
                               "Big Bed file (.bb) could not be created. Check if you have the bedToBigBed script in $PATH."],
            "MM_OUT_FOLDER_CREATION": [10, 0, "Could not create output folder."],
            "MM_NO_ARGUMENT": [11, 0, "Could not read the arguments. Make sure you provided an experimental matrix."],
            "MM_WRONG_EXPMAT": [12, 0,
                                "The experimental matrix could not be loaded. Check if it is correctly formatted and that your python version is >= 2.7."],
            "MM_WRONG_RANDPROP": [13, 0, "Proportion of random regions is too low (<= 0)."],
            "MM_LOW_NPROC": [14, 0, "Number of processors is too low (<= 0)."],
            "ME_OUT_FOLDER_CREATION": [15, 0, "Could not create output folder."],
            "ME_FEW_ARG": [16, 0,
                           "There are too few arguments. Please use -h option in order to verify command usage."],
            "ME_WRONG_EXPMAT": [17, 0,
                                "The experimental matrix could not be loaded. Check if it is correctly formatted and that your python version is >= 2.7."],
            "ME_MATCH_NOTFOUND": [18, 0,
                                  "Motif matching file was not found. Are you sure you performed --matching before --enrichment?"],
            "ME_BAD_MATCH": [19, 0, "Motif matching file is incorrect. Please perform --matching again."],
            "ME_LOW_NPROC": [20, 0, "Number of processor is too low (<= 0)."],
            "ME_RAND_NOTFOUND": [21, 0, "Random regions not found. Are you sure --matching was correctly performed?"],
            "ME_BAD_RAND": [22, 0, "Could not read random regions."],
            "ME_RAND_NOT_BED_BB": [23, 0, "Random regions are not in bed or bigbed format."],
            "MM_PSEUDOCOUNT_0": [24, 0,
                                 "There were errors involved in the creation of Position Weight Matrices. Some distributions of numpy  and/or scipy does not allow for pseudocounts == 0. Please increase the pseudocount (or use default value of 0.1) and try again."],
            "MM_MOTIFS_NOTFOUND": [25, 0, "A motif file was provided but it could not be loaded."],
            "XXXXXXX": [26, 0, "Xxxxxx"]
        }
        self.error_number = 0
        self.exit_status = 1
        self.error_message = 2

        self.warning_dictionary = {
            "DEFAULT_WARNING": [0, "Undefined warning."],
            "FP_ONE_REGION": [1,
                              "There are more than one 'regions' file in the experiment matrix. Only the first will be used."],
            "FP_MANY_DNASE": [2,
                              "There are more than one DNASE or ATAC 'reads' file. Only the first one will be used."],
            "FP_MANY_HISTONE": [3, "It is recomended that no more than three histone modifications should be used."],
            "FP_DNASE_PROC": [4, "The DNase  or ATAC file could not be processed."],
            "FP_HISTONE_PROC": [5, "The Histone file could not be processed."],
            "FP_SEQ_FORMAT": [6, "The DNase/ATAC+Histone sequence could not be formatted to be input for scikit."],
            "FP_HMM_APPLIC": [7, "The scikit HMM encountered errors when applied."],
            "MM_MANY_ARG": [8, "There are more than one arguments, only the first experimental matrix will be used."],
            "ME_MANY_ARG": [9, "There are more than two arguments, only the two first arguments will be used."],
            "ME_MANY_GENESETS": [10,
                                 "There are more than one geneset associated to a group, only the first one will be used."],
            "ME_FEW_GENESETS": [11, "There seems to be a geneset column, but no geneset was found for some entries."],
            "XXXXXXX": [12, "Xxxxxx"]
        }
        self.warning_number = 0
        self.warning_message = 1

    def throw_error(self, error_type, add_msg=""):
        """Throws the specified error type. If the error type does not exist, throws a default error message and exits.

        *Keyword arguments:*

            - error_type -- Error type.
            - add_msg -- Message to add to the error.
        """

        # Fetching error type
        try:
            error_number = self.error_dictionary[error_type][self.error_number]
            exit_status = self.error_dictionary[error_type][self.exit_status]
            error_message = self.error_dictionary[error_type][self.error_message]
        except (KeyError, IndexError):
            error_number = self.error_dictionary["DEFAULT_ERROR"][self.error_number]
            exit_status = self.error_dictionary["DEFAULT_ERROR"][self.exit_status]
            error_message = self.error_dictionary["DEFAULT_ERROR"][self.error_message]

        # Handling error
        complete_error_message = ("--------------------------------------------------\n"
                                  "Error Number: " + str(error_number) + ".\n"
                                                                         "Program: " + self.program_name + ".\n"
                                                                                                           "Report: " + error_message + " " + add_msg + "\n"
                                                                                                                                                        "Behaviour: The program will quit with exit status " + str(
            exit_status) + ".\n"
                           "--------------------------------------------------")
        print(complete_error_message, file=sys.stderr)
        traceback.print_exc()
        sys.exit(exit_status)

    def throw_warning(self, warning_type, add_msg=""):
        """Throws the specified warning type. If the warning type does not exist, throws a default warning message and exits.

        *Keyword arguments:*

            - warning_type -- Warning type.
            - add_msg -- Message to add to the error.
        """

        # Fetching warning type
        try:
            warning_number = self.warning_dictionary[warning_type][self.warning_number]
            warning_message = self.warning_dictionary[warning_type][self.warning_message]
        except (KeyError, IndexError):
            warning_number = self.warning_dictionary["DEFAULT_WARNING"][self.warning_number]
            warning_message = self.warning_dictionary["DEFAULT_WARNING"][self.warning_message]

        # Handling warning
        complete_warning_message = ("--------------------------------------------------\n"
                                    "Warning Number: " + str(warning_number) + ".\n"
                                                                               "Program: " + self.program_name + ".\n"
                                                                                                                 "Report: " + warning_message + " " + add_msg + "\n"
                                                                                                                                                                "--------------------------------------------------")
        print(complete_warning_message, file=sys.stderr)


###################################################################################################
# Auxiliary Functions as Static Methods
###################################################################################################

class AuxiliaryFunctions:
    """Class of auxiliary static functions."""

    @staticmethod
    def string_is_int(s):
        """Verifies if a string is a numeric integer.

        *Keyword arguments:*

            - s -- String to verify.
        """
        try:
            int(s)
            return True
        except ValueError:
            return False

    @staticmethod
    def string_is_float(s):
        """Verifies if a string is a numeric float.

        *Keyword arguments:*

            - s -- String to verify.
        """
        try:
            float(s)
            return True
        except ValueError:
            return False

    @staticmethod
    def correct_standard_bed_score(score):
        """Standardize scores between 0 and 1000.

        *Keyword arguments:*

            - score -- Score.
        """
        return min(max(score, 0), 1000)

    @staticmethod
    def overlap(t1, t2, strand_specific=False):
        """Checks if one interval contains any overlap with another interval.

        *Keyword arguments:*

            - t1 -- First tuple.
            - t2 -- Second tuple.
  
        *Return:*
            - -1 -- if i1 is before i2.
            - 1 -- if i1 is after i2.
            - 0 -- if there is any overlap.
        """
        if strand_specific:
            if t1[1] <= t2[0]: return -1  # interval1 is before interval2
            if t2[1] <= t1[0]: return 1  # interval1 is after interval2
            if t1[4] == t2[2]:
                return 0  # interval1 overlaps interval2
            else:
                return 2  # interval1 overlaps interval2 on the opposite strand
        else:
            if t1[1] <= t2[0]: return -1  # interval1 is before interval2
            if t2[1] <= t1[0]: return 1  # interval1 is after interval2
            return 0  # interval1 overlaps interval2

    @staticmethod
    def revcomp(s):
        """Revert complement string.

        *Keyword arguments:*

            - s -- String.
        """
        revDict = dict([("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"), ("N", "N")])
        return "".join([revDict[e] for e in s[::-1]])


"""
def strmatch(pattern, string, search="exact", case_insensitive=True):
    valid_types = ["exact", "inexact", "regex"]

    if case_insensitive:
        pattern = pattern.lower()
        string = string.lower()

    if search == "exact":
        return pattern == string
    elif search == "inexact":
        return pattern in string
    elif search == "regex":
        return re.search(pattern, string)
    else:
        raise ValueError("search must be one of these: {}".format(valid_types))


def cmp(a, b):
    return (a > b) - (a < b)


def npath(filename):
    ""Returns a normalised, absolute version of the path, with expanded user directory.""
    return os.path.abspath(os.path.expanduser(filename))


def get_rgtdata_path():
    return os.path.expanduser(os.getenv("RGTDATA", os.path.join(os.getenv("HOME"), "rgtdata")))


def which(program):
    ""Return path of program or None, see
    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python""

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None
    
"""

    
