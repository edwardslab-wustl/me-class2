# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 20:55:07 2018

@author: Chris
"""

#import numpy as np
class Sample:
    def __init__(self,sample_name,X,Y,L,E,Z,I):
        self.sample_name = sample_name
        self.individual_sample_names = set(sample_name.split("_"))
        ## Methylation Profiles        
        self.X = X 
        ## Expression Deciles           
        self.Y = Y 
        ## log10 gene length
        self.L = L 
        ## Number of exons
        self.E = E 
        ## 5hmC over TSS region
        self.Z = Z    
        ## Gene ID
        self.I = I

