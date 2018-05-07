# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 19:59:13 2018

@author: Chris
"""

import scipy
import numpy as np
import math


class Window:
    def __init__(self,low,high,bins,bp_window_size=None):
        self.low = low
        self.high = high
        self.len = self.high-self.low
        self.bins = bins
#        print "Window Bins: ",self.bins
        self.x = list()
        self.x_raw_val = list()
        self.y_val = list()
        self.y_raw_val = list()
        if bp_window_size == None:        
            self.set_bp_window_size()        
        else: 
            self.bp_window_size = bp_window_size   
#        print "BP Window Size: ",self.bp_window_size
            
    def clear_data(self):
        self.x = list()
        self.x_raw_val = list()
        self.y_val = list()
        self.y_raw_val = list()
            
    def set_counter(self,counter):
        self.counter = counter
    
    def set_bp_window_size(self):
        self.bp_window_size = (self.len)/self.bins

    def scale_sample_x_coords(self):
        self.scale_sample_x = [((x-self.low)/float(self.len))+self.counter for x in self.sample_x]
        
    def scale_raw_x_coords(self):
        self.scale_x_raw_val = [((x-self.low)/float(self.len))+self.counter for x in self.x_raw_val]

    def sample_window(self):
        if self.bp_window_size != 0:
            self.downsample_window()
        else:
            self.sample_x = list()
            self.sample_y_val = list()

    def downsample_window(self):
        self.sample_x = self.downsample(self.x,self.bp_window_size,self.bins)
#        print "Downsample: Number of bins in sample_x: ",len(self.sample_x)
        self.sample_y_val = self.downsample(self.y_val,self.bp_window_size,self.bins)
        
    def downsample(self,x,r,b):
        x = np.array(x)
        if (x.size % r) !=0:
            #if r < 3:
            #    r=3
            pad_size = math.ceil(float(x.size)/r)*r-x.size
            x_pad = np.append(x,np.zeros(pad_size)*np.NaN)
            #print("ps:", str(pad_size), "x_len:", str(x_pad.shape), "r:", str(r), "x:", str(x.size), file=sys.stderr)
            ds = scipy.nanmean(x_pad.reshape(-1,r),axis=1)
        else:
            ds = x.reshape(-1,int(r)).mean(axis=1)
            #ds = x.reshape(-1,r).mean(axis=1)
        #Average the remaining values if indivisible by bin size.
        if len(ds)>b:
            avg_lbs = np.mean(ds[b-1:])
            ds = ds[:b]
            ds[-1] = avg_lbs
        return list(ds)
