
def main():
    #
    ### Set up parser arguments
    global parser
    parser = setup_parser_arguments()
    ### Grab argument variables
    args = parser.parse_args()
    #
    dat_file_name = args.dat_file
    #
    with open(dat_file_name, 'r') as input_dat_file:
        input_dat_file_line = input_dat_file.readline()
    input_dat_file.close()
    #
    num_tss_wins = int( ( (input_dat_file_line.strip().split(','))[0] ).split(':')[1] )
    num_tss_bins = int( ( (input_dat_file_line.strip().split(','))[1] ).split(':')[1] )
    tss_offset   = int( ( (input_dat_file_line.strip().split(','))[2] ).split(':')[1] )
    bin_size = int(num_tss_wins / num_tss_bins)
    if args.verbose:
        sys.stderr.write("size of tss window: " + str( num_tss_wins ) + "\n" )
        sys.stderr.write("number of bins: " + str( num_tss_bins ) + "\n" )
        sys.stderr.write("tss_offset (0 means centered): " + str( tss_offset ) + "\n" )
        sys.stderr.write("bin size: " + str( bin_size ) + "\n" )
    #
    fi_file_lines = []
    for fi_file in args.file_names:
        if args.verbose:
            sys.stderr.write(fi_file + "\n")
        with open(fi_file, 'r') as input_file:
            fi_file_lines.append( input_file.readline() ) # reads only first line
        input_file.close()
    output_lines = []
    #
    total_num_bins = len( fi_file_lines[0].strip().split() )
    if total_num_bins % num_tss_bins != 0:
        sys.stderr.write("\nWARNING: total number of feature importances (" + str(total_num_bins) + ") is not divisible by the number of bins in the .dat file (" + str(num_tss_bins) + "). Will calculate averages anyway, but likely the .dat file is not from the same run as the featureImportance files. This will not affect the averages, but if the x-axis/bin locations could be wrong.\n\n")
    #
    for feat_shift in range( int( total_num_bins / num_tss_bins ) ):
        if feat_shift > 0:
            output_lines.append( "\n" )
        output_lines.append( "Feature" + str(feat_shift + 1) + "\n" )
        for tssbins in range( num_tss_bins ):
            #
            fi_value_ave = 0
            for line in fi_file_lines:
                fi_value_ave = fi_value_ave + float( (line.strip().split())[tssbins + feat_shift] )
            fi_value_ave = fi_value_ave / len(fi_file_lines)
            output_lines.append( str( int( (tssbins)*bin_size - (num_tss_wins/2 + tss_offset) ) ) + '\t' + str( fi_value_ave ) + '\n')

    # Write final output
#    with open('ave_feature_importance.txt', 'w') as output_file:
#        output_file.write(''.join(output_lines))
#    output_file.close()
    sys.stdout.write(''.join(output_lines))

def setup_parser_arguments():
    progDesc = '''This program takes a .dat file and corresponding featureImportance files use in an me-class2_classifier run and averages the feature importances across the different folds from cross-fold validation. Results are written to stdout as a tab-delimited file. The first column is the location relative to the tss of the feature. The second column is the average feature importance. If you have trouble, run with the --verbose flag for more information. If you use the 5mC_5hmC files, then Feature1 will be 5mC and Feature2 5hmC.'''
    progEpilog = '''Note: The .dat file is only used to specify the binning scheme for the x-axis values. If your .dat file is not being parsed correctly, you can create one by creating a file with the line:\n#tssWin:10000,tssBins:500,tssOffset:0\nReplace the window size, number of bins, and offset to match the values for your run. An offset of 0 corresponds to the TSS being centered in the window.  Save the file and specify this new file as the .dat file in the plotter.py command. If there are two features (e.g. 5mC and 5hmC) the file will have two sections, one for each file.\n\n\n'''
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,description=progDesc, epilog=progEpilog)
    ### Required positional arguments
    parser.add_argument('dat_file',help="one of the .dat files from the run, used to parse window and binning scheme")
    parser.add_argument('file_names',nargs='+',help="1 or more feature importance files to average together")
    parser.add_argument('-v','--verbose',default=False,action='store_true', \
            help="Verbose mode.  Prints file names and dat file info to stderr." )

    return parser

if __name__ == "__main__":
    import sys, os, argparse

    main()
