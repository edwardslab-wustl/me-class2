# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 16:42:29 2016

@author: cschlosberg and jredwards
"""

import meclass2.classifier
import sys

def main():
    ### Set up parser arguments
    #global parser
    parser = meclass2.classifier.setup_parser_arguments()
    ### Grab argument variables
    args = parser.parse_args()  
    samples = meclass2.classifier.load_samples(args.interp_dir,args)
    #sys.stderr.write("sample length: %s\n" % str(len(samples)))
    if len(samples) > 1:
        meclass2.classifier.loso_evaluation(samples,args)    
    else:
        sys.stderr.write("Only one sample, performing single sample cross validation\n")
        meclass2.classifier.crossfold_evaluation(samples, args)

if __name__ == '__main__':
    main()

