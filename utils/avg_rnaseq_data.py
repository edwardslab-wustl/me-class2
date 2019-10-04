
def main():
    #
    ### Set up parser arguments
    global parser
    parser = setup_parser_arguments()
    ### Grab argument variables
    args = parser.parse_args()
    #
    expr_info = {}
    for rnaseq_file in args.input_files:
        #
        if args.verbose:
            sys.stderr.write("reading in file: " + rnaseq_file + "\n")
        with open(rnaseq_file, 'r') as input_file:
            if args.skip_header:
                rnaseq_file_lines = input_file.readlines()[1:]
            else:
                rnaseq_file_lines = input_file.readlines()
        input_file.close()
        #
        if args.no_qual_check:
            qual_check = False 
        else:
            qual_check = True
        id_col = args.ID_col_number - 1
        expr_col = args.expr_col_number - 1
        for line in rnaseq_file_lines:
            line_items = line.strip().split()
            #if line.startswith('tracking_id') or (line_items[12] !='OK'):
            if qual_check and line_items[12] !='OK':
                continue
            try:
                expr_info[ line_items[id_col] ].append(float( line_items[expr_col] ))
            except KeyError:
                expr_info[ line_items[id_col] ] = []
                expr_info[ line_items[id_col] ].append(float( line_items[expr_col] ))
        #
        del rnaseq_file_lines[:]
         	
    output_lines = []
    if args.verbose:
        sys.stderr.write("Computing means...\n")
    for key in expr_info:
        output_lines.append( key + '\t' + str( mean( expr_info[key] ) ) + '\n')
        #output_lines.append( key + '\t' + str( np.mean( expr_info[key] ) ) + '\n')

    # Write final output
    if args.verbose:
        sys.stderr.write("Writing output...\n")
    if args.output:
        with open(args.output, 'w') as output_file:
            output_file.write(''.join(output_lines))
        output_file.close()
    else:
        sys.stdout.write(''.join(output_lines))


def mean( values ):
    total = 0
    number = 0
    for val in values:
        total += val
        number += 1
    mean = -1
    if number > 0:
        mean = total / number
	#sys.stderr.write("\tMEAN: " + str(mean) + "\n")
    else:
        if args.verbose:
            sys.stderr.write("\twarning, list has no values. returning a mean of -1\n")
    return mean

def setup_parser_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,description='''Computes the average expression amongst the provided files and outputs the GeneID and average expression. Output is to stdout unless otherwise specified.''')
    ### Required positional arguments
    parser.add_argument('input_files',nargs='+',help="1 or more RNA-seq results files")
    ### Optional arguments
    parser.add_argument('-o','--output',help="Output file. Otherwise results are printed to stdout.")
    parser.add_argument('-v','--verbose',action="store_true",help="Verbose output for troubleshooting.")
    parser.add_argument('--ID_col_number',default=4,help="Column that contains gene IDs. First column is 1. Default=4")
    parser.add_argument('--expr_col_number',default=10,help="Column that contains expression value. First column is 1. Default=10")
    parser.add_argument('--skip_header',action="store_true",help="Skip first line (header) in input files.")
    parser.add_argument('--no_qual_check',action="store_true",help='Do not perform qual check that column 13 value is "OK".')


    return parser

if __name__ == "__main__":
    import sys, argparse

    main()
