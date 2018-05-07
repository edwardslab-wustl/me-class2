# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 20:13:02 2018

@author: Chris
"""

class Sample:
    def __init__(self,sample_name, legend, predictions, expression, gene_ID):
        self.sample_name = sample_name
        #self.individual_sample_names = set(sample_name.split("_"))
        self.predictions = predictions
        self.expression = expression
        self.gene_ID = gene_ID    
        self.legend = legend

    def set_labels(self, labels):
        self.labels = labels