# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 20:09:30 2018

@author: Chris
"""

import sys

class FilterCounts:
    def __init__(self):
        ## Total number of entries        
        self.num_transcripts = 0
        ## Number of genes/transcripts defined as complete. cds start and end well defined.
        self.num_tx_complete = 0
        ## Number of genes/transcripts defined as incomplete. cds start and end ill defined.
        self.num_tx_incomplete = 0
        ## Number of unique gene ids
        self.num_uniq_gene_ids = 0 
        ## Number of genes ids with no transcripts after transcript filters
        self.num_gene_ids_no_tx_post_tx_filters = 0      
        ## Number of unique gene ids
        self.num_gene_ids_post_tx_filters = 0
        ## Number of genes on chromosomes defined in valid_chr    
        self.num_genes_simple_chrs = 0
        ## number of genes on chromosomes not defined in valid_chr    
        self.num_genes_other_chrs = 0
        ## Number of genes defined as protein coding (starting with NM_)
        self.num_genes_mrna = 0
         ## Number of genes defined as not protein coding (not starting with NM_)
        self.num_genes_notPC = 0
        ## Number of genes with unique TxStart sites
        self.num_genes_uniq_tx_start = 0
        ## Number of genes with unique TxStart sites
        self.num_genes_nonuniq_tx_start = 0
        ## Number of genes loaded into database for processing
        self.num_final_loaded_genes = 0
        
    def write_output(self,filter_flags):
        sys.stderr.write("Transcript Processing:\n")
        sys.stderr.write("\tIdentified %d transcripts\n" % (self.num_transcripts))
        if filter_flags.only_complete:
            sys.stderr.write("\tLoaded %d transcripts with well-defined CDS start and end\n" % (self.num_tx_complete))
            sys.stderr.write("\tLoaded %d transcripts with ill-defined CDS start and end\n" % (self.num_tx_incomplete))
        sys.stderr.write("Gene Processing:\n")        
        sys.stderr.write("\tIdentified %d unique gene annotations\n" % (self.num_uniq_gene_ids))
        if filter_flags.only_complete:
            sys.stderr.write("\tIdentified %d genes with no valid transcripts, post-transcript filter(s)\n" % (self.num_gene_ids_no_tx_post_tx_filters))
            sys.stderr.write("\tIdentified %d genes post-transcript filter(s)\n" % (self.num_gene_ids_post_tx_filters))
        if filter_flags.only_simple_chrs:
            sys.stderr.write("\tLoaded genes %d on simple chromosomes\n" % (self.num_genes_simple_chrs))        
            sys.stderr.write("\tLoaded genes %d on other chromosomes\n" % (self.num_genes_other_chrs))
        if filter_flags.only_mrna:
            sys.stderr.write("\tLoaded %d protein coding genes\n" % (self.num_genes_mrna))
            sys.stderr.write("\tLoaded %d non-protein coding genes\n" % (self.num_genes_notPC))
        if filter_flags.only_uniq_tx_start:
            sys.stderr.write("\tLoaded %d genes with unique Tx start site\n" % (self.num_genes_uniq_tx_start))
            sys.stderr.write("\tLoaded %d genes with non-unique Tx start site\n" % (self.num_genes_nonuniq_tx_start))
        sys.stderr.write("Loaded %d genes for cross-referencing\n" % (self.num_final_loaded_genes))
