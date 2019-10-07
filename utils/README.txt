README.txt

This is a set of utilities which may be useful for preprocessing data for me-class2 analysis.


######################
Processing bedms files
######################

There are three scripts provided to help manipulate .bedms files:
bedms_fix.pl  bedmspair2diff.pl  bedms_sane.pl

.bedms format is tab-delimited: chromsome, position1, position2, value
typically position2 = position1 + 1

bedms_sane.pl 
This script will check that a .bedms file is sorted. i.e. that every line for each chromosome is grouped together and within a chromosome all positions are in increasing order.

bedms_fix.pl
This script will sort a .bedms file by the criteria above.  Output is to stdout.

bedmspair2diff.pl
This script will take the difference between two .bedms files. If a position is missing in one file, it will be omitted from the output.



###################
Prep files for MLML
###################

fix_meth_files_for_mlml.py
MLML requires that each position with a value in the WGBS input file have a value in the TAB-seq (or oxBS-seq) file. This script will prune the input files to meet this criteria.


#################
Prep RNA-seq data
#################

avg_rnaseq_data.py
This script will take as input tab-delimited data files where at least one column is a gene ID and another is an expression column corresponding to multiple RNA-seq replicates and output a single file with the average methylation for each gene. This is a quick way to get values useful for fold-change calculations needed for ME-Class2.

combine_rnaseq.py
This script takes a refSeq annotation file (from the UCSC Genome Browser) and two expression files as input.  It converts the gene IDs into RefSeq IDs, floors the expression values, and collates expression data from the two files to create a properly formatted .expr expression data file for ME-Class2 analysis.


################################
Combine Feature Importance Files
################################

avg_featureImportance_files.py
This script takes multiple featureImportance files that are spit out using the --featureImportance flag of the meclass2_classifier, averages the values and computes the x-coordinate for each feature.

