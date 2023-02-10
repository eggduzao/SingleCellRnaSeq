#!/bin/bash

# Arguments
#input_matrix="--input_matrix=/home/egg/Projects/SingleCellRnaSeq/test/1_WAT/data/GSM4717152_CD45N-L1.csv"
#temporary_location="--temporary_location=/home/egg/Projects/SingleCellRnaSeq/test/1_WAT/temp/"
#output_location="--output_location=/home/egg/Projects/SingleCellRnaSeq/test/1_WAT/output/"

input_matrix="--input_matrix=/Users/egg/Projects/SingleCellRnaSeq/test/1_WAT/data/GSM4717152_CD45N-L1.csv"
temporary_location="--temporary_location=/Users/egg/Projects/SingleCellRnaSeq/test/1_WAT/temp/"
output_location="--output_location=/Users/egg/Projects/SingleCellRnaSeq/test/1_WAT/output/"

# Options
organism="--organism=hg38"
input_type="--input-type=csv"
cores="--cores=4"
output_type="--output-type=csv"
seed="--seed=123"

# Hidden Options
io_input_file_is_reversed="--io-input-file-is-reversed"


# Execution

scaw $organism $input_type $cores $output_type $seed $io_input_file_is_reversed $input_matrix $temporary_location $output_location 


