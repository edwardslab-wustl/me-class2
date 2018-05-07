# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 20:05:52 2018

@author: Chris
"""

class Transcript:
    def __init__(self,transcript_id,chrom,strand,pos_low,pos_high,exons,common_name,cds_tx_start_cmpl,cds_tx_end_cmpl):
        self.transcript_id = transcript_id
        self.chrom = chrom
        self.strand = strand
        self.pos_low = pos_low
        self.pos_high = pos_high
        self.exons = exons
        self.common_name = common_name
        self.cds_tx_start_cmpl = cds_tx_start_cmpl
        self.cds_tx_end_cmpl = cds_tx_end_cmpl