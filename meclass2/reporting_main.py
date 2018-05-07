# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 20:23:54 2018

@author: Chris and John
"""
#Next two lines are so does not use X-windows for display.
import matplotlib as mpl 
mpl.use('Agg')
import meclass2.reporting
import copy


def main():
    ### Set up parser arguments
    #global parser
    parser = meclass2.reporting.setup_parser_arguments()
    ### Grab argument variables
    args = parser.parse_args()  
    expression_dict,gene_ID_dict,gene_name_dict,samples_dict = meclass2.reporting.load_labels(args.pred_files, args)
    samples = meclass2.reporting.load_samples(args.pred_files,expression_dict,gene_ID_dict,samples_dict,args)
    meclass2.reporting.calculate_labels(samples, args)
    results = meclass2.reporting.calculate_results(samples,expression_dict,gene_ID_dict,gene_name_dict,samples_dict,args)
    results_P = copy.deepcopy(results)
    results_N = copy.deepcopy(results)
    meclass2.reporting.calculate_roc(results,args)
    meclass2.reporting.print_roc(results, args)
    meclass2.reporting.calculate_acc_rejectrate(results,'all',args)
    meclass2.reporting.calculate_acc_rejectrate(results_P,'pos',args)
    meclass2.reporting.calculate_acc_rejectrate(results_N,'neg',args)
    meclass2.reporting.print_acc_rejectrate(results, 'all', args)
    meclass2.reporting.print_acc_rejectrate(results_P, 'pos', args)
    meclass2.reporting.print_acc_rejectrate(results_N, 'neg', args)
    meclass2.reporting.print_high_acc_genes(results, args)
    if args.plot_results:
        version = int(mpl.__version__.split('.')[0])
        meclass2.reporting.set_matplotlib_params(args, version)
        outFile = args.outFileBase + ".roc.png"
        meclass2.reporting.plot_roc(results, outFile, args, version)
        outFile = args.outFileBase + ".all.acc_rejectrate.png"
        meclass2.reporting.plot_acc_rejectrate(results, outFile, args, version)
        outFile = args.outFileBase + ".pos.acc_rejectrate.png"
        meclass2.reporting.plot_acc_rejectrate(results_P, outFile, args, version)
        outFile = args.outFileBase + ".neg.acc_rejectrate.png"
        meclass2.reporting.plot_acc_rejectrate(results_N, outFile, args, version)


if __name__ == '__main__':
    main()

