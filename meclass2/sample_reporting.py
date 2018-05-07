# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 20:55:07 2018

@author: Chris
"""

class Sample:
    def __init__(self,sample_name, legend, predictions, expression, gene_ID,
            sample_names_in_label_file):
        self.sample_name = sample_name
        self.sample_names_in_label_file = sample_names_in_label_file
        #self.individual_sample_names = set(sample_name.split("_"))
        self.predictions = predictions
        self.expression = expression
        self.gene_ID = gene_ID    
        self.legend = legend

    def set_labels(self, labels):
        self.labels = labels

