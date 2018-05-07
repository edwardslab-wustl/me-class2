# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 20:34:25 2018

@author: Chris
"""

import sys
#Next two lines are so does not use X-windows for display when plotting
import matplotlib as mpl 
mpl.use('Agg')
import meclass2.interpolation

### For help call: python <script>.py -h

def main():

    ### Set up parser arguments
    #global parser
    parser = meclass2.interpolation.setup_parser_arguments()
    ### Grab argument variables
    args = parser.parse_args()
    if args.interp_type == "HRPS":
        args.interp_type = "TSS"
    ### Define valid chromosomes
    global valid_chr
    if args.autosomeFlag:
        valid_chr = set(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8',
            'chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16',
            'chr17','chr18','chr19','chr20','chr21','chr22'])
    else:
        valid_chr = set(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8',
            'chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16',
            'chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'])
    ### Load gene annotations
    gene_annot,filter_flags = meclass2.interpolation.load_gene_annot(args.annot_file,valid_chr)
    ### Load sample files 
    sample_files = meclass2.interpolation.load_sample_info(args)    
    ### Load gene expression
    gene_list,differential = meclass2.interpolation.load_gene_expression_rpkm(sample_files,gene_annot)
    ### Release memory from holding ENSG database
    del gene_annot      
    ### Print Label List
    label_list_fh = meclass2.interpolation.print_label_list_header(gene_list,args.tag) 
    if args.interp_cpg_density:    
        ### Load up all CpG positions
        cpg_density = dict()
        cpg_density = meclass2.interpolation.load_cpg_perc_ref(args.cpg_density_file)   
    else:
        cpg_density = None
    if args.interp_type=="DMR":
        meclass2.interpolation.process_dmr(sample_files,gene_list,label_list_fh,differential,args)
    else:
        ### Calculate and print interpolation for each gene
        meclass2.interpolation.interpolate_methylation_cpgdensity_windows(sample_files,gene_list,cpg_density,label_list_fh,differential,valid_chr,args)    
    sys.stderr.write("Done!\n")

if __name__ == '__main__':
    main()
