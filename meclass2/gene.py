# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 20:00:21 2018

@author: Chris
"""

class Gene:
    def __init__(self,gene_id):
        ### Either the REFSEQ or ENSG ID        
        self.gene_id = gene_id
        ### List of all the transcripts
        self.transcripts = list()
        
    def set_location_information(self,chrom,strand,pos_low,pos_high,exons,common_name):        
        self.chrom = chrom
        self.strand = strand
        self.pos_low = pos_low
        self.pos_high = pos_high
        self.exons = exons
        self.common_name = common_name

    def set_half_tss_window(self,half_tss_window):
        self.half_tss_window = half_tss_window

    def set_expr(self,expr):
        self.expr = expr    
    
    def set_sample(self,sample):
        self.sample = sample
        
    def append_transcript(self,transcript):
        self.transcripts.append(transcript)

    def set_introns(self):
        self.introns = self.calculate_introns()        
        
    def calculate_introns(self):
        introns = list()    
        for i,exon_tup in enumerate(self.exons):
            try:
                i_tup = (exon_tup[1]+1,self.exons[i+1][0]-1) 
                introns.append(i_tup)
            except IndexError:
                continue
        return introns        
        
    def merge_gene_transcripts(self):
        common_names = list()        
        for i,t in enumerate(self.transcripts):
            common_names.append(t.common_name)            
            if i == 0:
                merged_exons = merge_intervals(t.exons)
            else:
                current_merged_exons = merged_exons+t.exons
                merged_exons = merge_intervals(current_merged_exons)                              
        pos_low = merged_exons[0][0]
        pos_high = merged_exons[-1][1]
        chrom = self.transcripts[0].chrom
        strand = self.transcripts[0].strand
        common_name = get_most_common_element(common_names)
        self.set_location_information(chrom,strand,pos_low,pos_high,merged_exons,common_name)
        self.set_introns()
        return