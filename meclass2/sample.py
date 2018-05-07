# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 20:55:07 2018

@author: Chris
"""

import numpy as np

class Sample:
    def __init__(self,sample_name,X,Y,C,P,M,E,Z,geneIDs,data_info):
        self.sample_name = sample_name
        self.X = np.asarray(X)        
        self.Y = np.asarray(Y)
        self.C = np.asarray(C)
        self.P = P  
        self.M = np.asarray(M)
        self.E = np.asarray(E)
        self.Z = np.asarray(Z)
        self.geneIDs = geneIDs
        self.data_info = data_info # tssWin, tssBins, tssOffset)
        
