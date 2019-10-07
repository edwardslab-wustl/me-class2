
def main():
    #
    ### Set up parser arguments
    global parser
    parser = setup_parser_arguments()
    ### Grab argument variables
    args = parser.parse_args()
    #
    df_gene_annot = pd.read_csv(args.refGene_file, sep='\t', header=None, comment='#', na_values='NA', \
                    names = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames'])
    #annot_columns = df_gene_annot.columns
    #
    # Drop which are cdsStartStat != cmpl; cdsEndStat != cmpl
    index_to_drop = df_gene_annot[ (df_gene_annot['cdsStartStat'] != 'cmpl') |  (df_gene_annot['cdsEndStat'] != 'cmpl') ].index
    df_gene_annot.drop(index_to_drop , inplace=True)

    #Remove duplicate gene ids by keeping lowest accession
    #
    df_gene_annot = df_gene_annot.set_index('name2', drop=False)
    df_gene_annot[['forsort1', 'forsort2']] = df_gene_annot.name.str.split('_', expand=True)
    df_gene_annot['forsort2'] = df_gene_annot['forsort2'].astype(int)
    df_gene_annot = df_gene_annot.sort_values(by=['forsort2'])
    df_gene_annot = df_gene_annot.drop_duplicates(subset='name2', keep='first')
    #
    # Readin expression files
    #
    df1_expr = pd.read_table(args.file1, index_col=False, header=None, na_values='NA', names = ['gene_id_1', 'ave_expr_1'])
    df1_expr = df2_expr.set_index('gene_id_1')
    #
    df2_expr = pd.read_table(args.file2, index_col=False, header=None, na_values='NA', names = ['gene_id_2', 'ave_expr_2'])
    df2_expr = df1_expr.set_index('gene_id_2')
    #
    if not args.turn_off_floor:
        #df1_expr['ave_expr_1']=df1_expr['ave_expr_1'].apply(floor)
        #df2_expr['ave_expr_2']=df2_expr['ave_expr_2'].apply(floor)
        df1_expr['ave_expr_1']=df1_expr['ave_expr_1'].apply(floor)
        df2_expr['ave_expr_2']=df2_expr['ave_expr_2'].apply(floor)

        # Alternate way
        # df1_expr.loc[df1_expr['ave_expr_1'] < 5, 'ave_expr_1'] = 5
        # df2_expr.loc[df2_expr['ave_expr_2'] < 5, 'ave_expr_2'] = 5
        #OR *Better*
        # df1_expr['ave_expr_1'] = np.where( df1_expr['ave_expr_1'] < 5, 5, df1_expr['ave_expr_1'] )	
        # df2_expr['ave_expr_2'] = np.where( df2_expr['ave_expr_2'] < 5, 5, df2_expr['ave_expr_2'] )

	
    df_merged = pd.concat( [df_gene_annot, df1_expr, df2_expr], axis=1, join='inner') 
    #
    # Write output
    df_merged[['name', 'ave_expr_1', 'ave_expr_2']].to_csv(args.output, sep='\t', na_rep='NA', header=None, index=False)


def floor( value ):
    floor = 5
    if value < floor:
        new_val = floor
    else: 
        new_val = value
    return new_val

def setup_parser_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,description='''Converts gene symbols to RefSeq IDs and collates expression data from two samples.''')
    ### Required positional arguments
    parser.add_argument('refGene_file',help="Refseq annotation file from UCSC genome browser.")
    parser.add_argument('file1',help="RNA-seq results file from first sample.")
    parser.add_argument('file2',help="RNA-seq results file from second sample.")
    parser.add_argument('output',help="Name of output file.")
    ### Optional arguments
    parser.add_argument('-o','--output',help="Output file. Otherwise results are printed to stdout.")
    parser.add_argument('-v','--verbose',action="store_true",help="Verbose output for troubleshooting.")
    parser.add_argument('--floor',default=5,help="Set expression values lower than this to this value. Default=5")
    parser.add_argument('--turn_off_floor',action="store_true",help="Do not floor the expression values.")



    return parser

if __name__ == "__main__":
    import sys, argparse
    import pandas as pd
    import numpy as np

    main()




