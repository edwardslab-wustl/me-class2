# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 20:07:59 2018

@author: Chris
"""

class FilterFlags:
    def __init__(self,annot_filetype="",only_simple_chrs=True,only_mrna=True,only_complete=True,only_uniq_tx_start=False,differential=False,collapse_on_common_name=True,exclude_lt_four_exons=False):
        ### Gene Annotation Format 
        self.annot_filetype = annot_filetype
        ### Transcript Level Filters        
        ## Only look at transcripts defined as complete. cds start and end well defined.
        self.only_complete = only_complete
        ### Gene Level Filters
        ## Only look at genes on chromosomes defined in valid_chr    
        self.only_simple_chrs = only_simple_chrs
        ## Only look at genes defined as protein coding (starting with NM_)
        self.only_mrna = only_mrna
        ## Only look at genes with unique TxStart sites
        self.only_uniq_tx_start = only_uniq_tx_start
        ## Identify that this interpolation is for differential
        self.differential = differential
        ## Collapse Annotations so that only common gene name represented
        self.collapse_on_common_name = collapse_on_common_name
        ## Exclude genes with < 4 exons to have straightforward comparison with ROI
        self.exclude_lt_four_exons=exclude_lt_four_exons
