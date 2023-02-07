#!/bin/bash

# Arguments
input_matrix="--input_matrix=/home/egg/Projects/SingleCellRnaSeq/test/2_CellRangerInput/raw_feature_bc_matrix/"
temporary_location="--temporary_location=/home/egg/Projects/SingleCellRnaSeq/test/2_CellRangerInput/temp/"
output_location="--output_location=/home/egg/Projects/SingleCellRnaSeq/test/2_CellRangerInput/output/"

#input_matrix="--input_matrix=/Users/egg/Projects/SingleCellRnaSeq/test/data/GSM4717152_CD45N-L1.csv"
#temporary_location="--temporary_location=/Users/egg/Projects/SingleCellRnaSeq/test/temp/"
#output_location="--output_location=/Users/egg/Projects/SingleCellRnaSeq/test/output/"

# Options
organism="--organism=hg38"
input_type="--input-type=cellranger"
cores="--cores=4"
output_type="--output-type=csv"
seed="--seed=123"

# Hidden Options
io_input_file_is_reversed="--io-input-file-is-reversed" # $io_input_file_is_reversed


# Execution

scaw $organism $input_type $cores $output_type $seed $input_matrix $temporary_location $output_location 


