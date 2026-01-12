#!/usr/bin/env Rscript
#Given one list of samples, split [1] in two files with the samples indicated in [2]

#[1] First argument: input file that we want to split
#[2] Second argument: list of samples of the first condition
#[3] Third argument: output file of the first condition

# Parse command line arguments
print("Parsing samples...")
CHARACTER_command_args <- commandArgs(trailingOnly=TRUE)

# Load the input file
print(paste0("Loading ",CHARACTER_command_args[1],"..."))
input_file <- read.table(CHARACTER_command_args[1],header=TRUE)

# Load the list of samples of the first condition
# Replace the dashes with dots
formatted_string1 <- gsub("-",".",CHARACTER_command_args[2])
first_condition <- unlist(strsplit(formatted_string1,","))

# Take the samples of first condition and generate a file with just these columns
stopifnot(first_condition %in% colnames(input_file))
first_output <- input_file[first_condition]

# Save the output files
path1 <- CHARACTER_command_args[3]

print(paste0("Writing ",path1))
write.table(first_output,file=path1,quote=FALSE,sep="\t")
print(paste0("Saved ",path1))
