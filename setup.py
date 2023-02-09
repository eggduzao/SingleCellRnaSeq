from __future__ import print_function
"""
Setup Module
===================
Placeholder.

Authors: Eduardo G. Gusmao. 

Install Command: pip install --user .
                 pip install --user . --upgrade
"""

# Python
import os
import io
import re
import sys
import pwd
import shutil
import optparse
import setuptools

###################################################################################################
# Parameters
###################################################################################################

"""
Tools Dictionary:
  * Insert in the dictionary bellow a key+tuple X: (Y,Z,W,K) representing:
    - X: A string representing the name the user should provide for this script in order to install the tool.
    - Y: A string representing the name of the program which the user should type in the terminal in order to execute the tool.
    - Z: A string representing the path to the main function that executes the program.
    - W: A list with the package requirements for that program.
    - K: A list with external binary files that will be copied by the setup installation function to the user's bin folder.
Tools Dictionary Standard:
  * All programs should start with "bk-" followed by the program name.
  * The main function called within the script must be termed "main".
"""

# Common dependencies.
common_deps = ["numpy>=1.21.0",
               "scipy>=1.10.0"]

tools_dictionary = {
    "scaw": (
        "scaw",
        "src.main:main",
        ["numpy>=1.21.0", "scipy>=1.10.0", "pandas>=1.5.0", "scanpy>=1.9.0", "anndata>=0.8.0"],
        []
    )
}


###################################################################################################
# Unsupported Platforms
###################################################################################################

supported_platforms = ["linux", "linux2", "darwin"]
if sys.platform not in supported_platforms:
    print("ERROR: This package currently supports only unix-based systems (Linux).")
    exit(0)


###################################################################################################
# Auxiliary Functions/Classes
###################################################################################################

# PassThroughOptionParser Class
class PassThroughOptionParser(optparse.OptionParser):
    """
    An 'unknown option' pass-through implementation of OptionParser.
    When unknown arguments are encountered, bundle with largs and try again,
    until rargs is depleted.
    sys.exit(status) will still be called if a known argument is passed
    incorrectly (e.g. missing arguments or bad argument types, etc.)
    """

    def _process_args(self, largs, rargs, values):
      while rargs:
        try:
          optparse.OptionParser._process_args(self, largs, rargs, values)
        except (optparse.BadOptionError, optparse.AmbiguousOptionError) as err:
          largs.append(err.opt_str)

# recursive_chown_chmod Function
def recursive_chown_chmod(path_to_walk, uid, gid, file_permission, path_permission):
    """
    Recursively applies chown from path.
    """
    for root_dir, directory_list, file_list in os.walk(path_to_walk):
        os.chown(root_dir, uid, gid)
        os.chmod(root_dir, path_permission)
        for f in file_list:
            current_complete_file = os.path.join(root_dir, f)
            os.chown(current_complete_file, uid, gid)
            os.chmod(current_complete_file, file_permission)

def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")
  
  
###################################################################################################
# Processing Input Arguments
###################################################################################################

# Parameters
current_version = find_version("src", "__version__.py")
usage_message = "python setup.py install [python options] [SCAW options]"
version_message = "SCAW. Version: " + str(current_version)

# Initializing Option Parser
parser = PassThroughOptionParser(usage = usage_message, version = version_message)

# Parameter: Copy scaw data folder
param_copy_scaw_data_name = "--copy-scaw-data"
param_copy_scaw_data_help = "Explain here."
parser.add_option(param_copy_scaw_data_name, dest="param_copy_scaw_data", action="store_false", default=True, help=param_copy_scaw_data_help)

# Processing Options
options, arguments = parser.parse_args()
param_copy_scaw_data = options.param_copy_scaw_data

# Manually Removing Additional Options from sys.argv
new_sys_argv = []
for e in sys.argv:
    if param_copy_scaw_data_name == e[:len(param_copy_scaw_data_name)]:
        continue
    new_sys_argv.append(e)
sys.argv = new_sys_argv

# Defining entry points
current_entry_points = {"console_scripts": []}
for tool_option in tools_dictionary.keys():
    current_entry_points["console_scripts"].append(" = ".join(tools_dictionary[tool_option][:2]))

# Defining install requirements
current_install_requires = common_deps
for tool_option in tools_dictionary.keys():
    current_install_requires += tools_dictionary[tool_option][2]


###################################################################################################
# Creating Data Path
###################################################################################################

# Default scaw_data_path
current_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
scaw_data_path = os.path.join(os.path.expanduser("~"), "scaw_data")

# Creating data path
if(param_copy_scaw_data):
    if(os.path.exists(scaw_data_path)):
        shutil.rmtree(scaw_data_path)
    shutil.copytree(current_path, scaw_data_path)

# Creating data.config
data_config_file_name = os.path.join(scaw_data_path, "data.config")
# if not os.path.isfile(data_config_file_name):
data_config_file = open(data_config_file_name, "w")
data_config_file.write("# Configuration file loaded at SCAW startup. CAREFUL: any changes shall be overwritten\n"
                       "# whenever SCAW is (re)installed.\n\n")

genome = "mm9"
data_config_file.write("[" + genome + "]\n")
data_config_file.write("genome: " + os.path.join(genome, "genome_mm9.fa\n"))
data_config_file.write("chromosome_sizes: " + os.path.join(genome, "chrom.sizes.mm9\n"))
data_config_file.write("gene_regions: " + os.path.join(genome, "genes_Gencode_mm9.bed\n"))
data_config_file.write("# gene_regions: " + os.path.join(genome, "genes_RefSeq_mm9.bed # alternative to Gencode\n"))
data_config_file.write("annotation: " + os.path.join(genome, "gencode.vM1.annotation.gtf\n"))
data_config_file.write("gene_alias: " + os.path.join(genome, "alias_mouse.txt\n\n"))
data_config_file.write("repeat_maskers: " + os.path.join(genome, "repeat_maskers\n\n"))

genome = "mm10"
data_config_file.write("[" + genome + "]\n")
data_config_file.write("genome: " + os.path.join(genome, "genome_mm10.fa\n"))
data_config_file.write("chromosome_sizes: " + os.path.join(genome, "chrom.sizes.mm10\n"))
data_config_file.write("gene_regions: " + os.path.join(genome, "genes_Gencode_mm10.bed\n"))
data_config_file.write("# gene_regions: " + os.path.join(genome, "genes_RefSeq_mm10.bed # alternative to Gencode\n"))
data_config_file.write("annotation: " + os.path.join(genome, "gencode.vM25.annotation.gtf\n"))
data_config_file.write("gene_alias: " + os.path.join(genome, "alias_mouse.txt\n\n"))

genome = "mm39"
data_config_file.write("[" + genome + "]\n")
data_config_file.write("genome: " + os.path.join(genome, "genome_mm39.fa\n"))
data_config_file.write("chromosome_sizes: " + os.path.join(genome, "chrom.sizes.mm39\n"))
data_config_file.write("gene_regions: " + os.path.join(genome, "genes_Gencode_mm39.bed\n"))
data_config_file.write("annotation: " + os.path.join(genome, "gencode.vM25.annotation.gtf\n"))
data_config_file.write("gene_alias: " + os.path.join(genome, "alias_mouse.txt\n\n"))

genome = "hg19"
data_config_file.write("[" + genome + "]\n")
data_config_file.write("genome: " + os.path.join(genome, "genome_hg19.fa\n"))
data_config_file.write("chromosome_sizes: " + os.path.join(genome, "chrom.sizes.hg19\n"))
data_config_file.write("gene_regions: " + os.path.join(genome, "genes_Gencode_hg19.bed\n"))
data_config_file.write("# gene_regions: " + os.path.join(genome, "genes_RefSeq_hg19.bed # alternative to Gencode\n"))
data_config_file.write("annotation: " + os.path.join(genome, "gencode.v19.annotation.gtf\n"))
data_config_file.write("gene_alias: " + os.path.join(genome, "alias_human.txt\n\n"))
data_config_file.write("repeat_maskers: " + os.path.join(genome, "repeat_maskers\n\n"))

genome = "hg38"
data_config_file.write("[" + genome + "]\n")
data_config_file.write("genome: " + os.path.join(genome, "genome_hg38.fa\n"))
data_config_file.write("chromosome_sizes: " + os.path.join(genome, "chrom.sizes.hg38\n"))
data_config_file.write("gene_regions: " + os.path.join(genome, "genes_Gencode_hg38.bed\n"))
data_config_file.write("# gene_regions: " + os.path.join(genome, "genes_RefSeq_hg38.bed # alternative to Gencode\n"))
data_config_file.write("annotation: " + os.path.join(genome, "gencode.v21.annotation.gtf\n"))
data_config_file.write("gene_alias: " + os.path.join(genome, "alias_human.txt\n\n"))
data_config_file.write("repeat_maskers: " + os.path.join(genome, "repeat_maskers\n\n"))

genome = "zv9"
data_config_file.write("[" + genome + "]\n")
data_config_file.write("genome: " + os.path.join(genome, "genome_zv9_ensembl_release_79.fa\n"))
data_config_file.write("chromosome_sizes: " + os.path.join(genome, "chrom.sizes.zv9\n"))
data_config_file.write("gene_regions: " + os.path.join(genome, "genes_zv9.bed\n"))
data_config_file.write("annotation: " + os.path.join(genome, "Danio_rerio.Zv9.79.gtf\n"))
data_config_file.write("gene_alias: " + os.path.join(genome, "alias_zebrafish.txt\n\n"))

######

session = "Cooler"
data_config_file.write("[" + session + "]\n")
data_config_file.write("command: cooler\n\n")

session = "Juicer"
data_config_file.write("[" + session + "]\n")
data_config_file.write("command: java\n")
data_config_file.write("options: -Djava.awt.headless=true -Xmx32000m -jar\n")
data_config_file.write("jar: bin/juicer_tools_1.22.01.jar\n\n")

session = "Report"
data_config_file.write("[" + session + "]\n")
data_config_file.write("texsettingsandmaterials: report/texsettingsandmaterials\n")
data_config_file.write("reportstylebasic: report/reportstylebasic.sty\n")
data_config_file.write("reportstylefull: report/reportstylefull.sty\n")
data_config_file.write("executable: report/run.sh\n\n")


data_config_file.close()


###################################################################################################
# Setup Function
###################################################################################################

# Parameters
package_name = "SCAW"
package_version = str(current_version)
short_description = "Computational Workflow to analyze single-cell RNA-seq data"
short_description = "Toolkit to perform regulatory genomics data analysis"
classifiers_list = ["Topic :: Scientific/Engineering :: Bio-Informatics",
                    "Topic :: Scientific/Engineering :: Artificial Intelligence"]
keywords_list = ["scRNA-seq"]
author_list = ["Eduardo G. Gusmao"]
corresponding_mail = "eduardo.gusmao@utu.fi"
license_type = "GPL"

# Fetching additional structural files
readme_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "README.md")

# Fetching long description
readme_file = open(readme_file_name, "rU", encoding="utf-8")
long_description = readme_file.read()
readme_file.close()

# Setup Function
setuptools.setup(name = package_name,
                 version = package_version,
                 description = short_description,
                 long_description = long_description,
                 classifiers = classifiers_list,
                 keywords = ", ".join(keywords_list),
                 author = ", ".join(author_list),
                 author_email = corresponding_mail,
                 license = license_type,
                 packages = setuptools.find_packages(),
                 entry_points = current_entry_points,
                 install_requires = current_install_requires
)


###################################################################################################
# Termination
###################################################################################################

# Modifying Permissions when Running Superuser/Admin
if(param_copy_scaw_data):

    # Get current user and set default permissions for files to be visible and binaries executable
    current_user = os.getenv("SUDO_USER")
    default_file_permission = 0o644
    default_path_permission = 0o755

    # Setting the permissions
    if current_user:
        current_user_uid = pwd.getpwnam(current_user).pw_uid
        current_user_gid = pwd.getpwnam(current_user).pw_gid
        recursive_chown_chmod(scaw_data_path, current_user_uid, current_user_gid, default_file_permission,
                              default_path_permission)



