

def main(): 
    parser = setup_parser_arguments()
    ### Grab argument variables
    args = parser.parse_args()

    file1 = args.file1
    file2 = args.file2
    if args.verbose:
        sys.stderr.write("File1: " + file1 + "\n")
        sys.stderr.write("File2: " + file2 + "\n")

    df1_bed = pd.read_table(file1, index_col=False, header=None, skiprows=1, na_values='NA', names = ['chrom1', 'start1', 'strand1', 'context1', 'ratio1', 'tot_read1'])
    df1_bed.sort_values(by=['chrom1', 'start1']) # sort (doesn't change answer)
    df1_bed['idx_value1'] = df1_bed['chrom1']+'_'+df1_bed['start1'].astype(str)
    df1_bed = df1_bed.set_index('idx_value1')

    ########

    df2_bed = pd.read_table(file2, index_col=False, header=None, skiprows=1, na_values='NA', names = ['chrom2', 'start2', 'strand2', 'context2', 'ratio2', 'tot_read2'])
    df2_bed.sort_values(by=['chrom2', 'start2']) # sort (doesn't change answer)
    df2_bed['idx_value2'] = df2_bed['chrom2']+'_'+df2_bed['start2'].astype(str)
    df2_bed = df2_bed.set_index('idx_value2')
    #
    df_merged = df_merged = pd.concat( [ df1_bed, df2_bed ], axis=1, join='inner')
    #
    cols_to_keep_1 = ['chrom1', 'start1', 'strand1', 'context1', 'ratio1', 'tot_read1']
    cols_to_keep_2 = ['chrom2', 'start2', 'strand2', 'context2', 'ratio2', 'tot_read2']	
    # write output
    df_merged[cols_to_keep_1].to_csv(file1 + '.fixed', sep='\t', na_rep='NA', header=None, index=False)
    df_merged[cols_to_keep_2].to_csv(file2 + '.fixed', sep='\t', na_rep='NA', header=None, index=False)
    #---------------------

def setup_parser_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,description='''
Reformat methylation data for MLML. Input data should be tab-delimited 4 column: Chomosome<tab>position<tab>meth_level<tab>Coverage. This program will prune both supplied files so only positions with data in each file are output, as required by MLML.''')
    ### Required positional arguments
    parser.add_argument('file1',default="",help="First .meth file.")
    parser.add_argument('file2',default="",help="Second .meth file.")
    parser.add_argument('-v','--verbose',action='store_true',help="Verbose output for troubleshooting.")

    return parser

if __name__ == "__main__":
    import sys, os, argparse
    import pandas as pd
    main()
