# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 21:02:29 2018

@author: Chris
"""

#Next two lines are so does not use X-windows for display.
import matplotlib as mpl 
mpl.use('Agg') 
import meclass2.clustering

def main():
    ### Set up parser arguments
    global parser
    parser = meclass2.clustering.setup_parser_arguments()
    ### Grab argument variables
    args = parser.parse_args()
    #low_pred_bound = 0.7
    #low_pred_bound = 0.8
    #high_pred_bound = 1.0
    low_pred_bound = float(args.lowerPredBound)
    high_pred_bound = float(args.upperPredBound)
    corr_pred = 1
    samples = meclass2.clustering.process_samples(args.interp_dir,args.tag,args) 
    meclass2.clustering.all_sample_clustering(args.base,samples,low_pred_bound,high_pred_bound,corr_pred,int(args.upperBound),int(args.lowerBound),int(args.numClusters),args)
 
if __name__ == '__main__':
    main()
