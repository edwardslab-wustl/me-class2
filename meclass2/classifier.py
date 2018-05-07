# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 16:42:29 2016

@author: cschlosberg
"""

### Custom container Class definitions
from meclass2.sample_classifier import Sample
import sys, os, argparse, linecache, operator, math, pickle, collections, random, copy
import re
import numpy as np
import multiprocessing
import scipy.stats
import sklearn.model_selection
import sklearn.ensemble
import sklearn.metrics
import sklearn.neighbors    
import sklearn.externals


##  class Sample:
##      def __init__(self,sample_name,X,Y,L,E,Z,I):
##          self.sample_name = sample_name
##          self.individual_sample_names = set(sample_name.split("_"))
##          ## Methylation Profiles        
##          self.X = X 
##          ## Expression Deciles           
##          self.Y = Y
##          ## log10 gene length
##          self.L = L
##          ## Number of exons
##          self.E = E
##          ## 5hmC over TSS region
##          self.Z = Z    
##          ## Gene ID
##          self.I = I         

## def main():
##     ### Set up parser arguments
##     global parser
##     parser = setup_parser_arguments()
##     ### Grab argument variables
##     args = parser.parse_args()  
##     samples = load_samples(args.interp_dir)
##     #sys.stderr.write("sample length: %s\n" % str(len(samples)))
##     if len(samples) > 1:
##         loso_evaluation(samples,args)    
##     else:
##         sys.stderr.write("Only one sample, performing single sample cross validation\n")
##         crossfold_evaluation(samples, args)



### Function Definitions below here

def crossfold_evaluation(samples, args):
    X_data = list()
    Y_data = list()
    Z_data = list()
    I_data = list()
    sample_name = ''
    for i,sample in enumerate(samples):
        sys.stderr.write("Preparing evaluation for %s\n" % (sample.sample_name))
        X_data+=sample.X
        Y_data+=sample.Y
        Z_data+=sample.Z
        I_data+=sample.I
        sample_name=sample.sample_name
    total_genes = len(I_data)
    X_data = np.asarray(X_data)
    Y_data = np.asarray(Y_data)
    Z_data = np.asarray(Z_data)
    I_data = np.asarray(I_data)
    sys.stderr.write("total genes: %s\n" % (total_genes))
    n_folds = int(args.k_folds)
    #method = sklearn.ensemble.RandomForestClassifier(n_estimators=2001,n_jobs=args.num_jobs)    
    method = sklearn.ensemble.RandomForestClassifier(n_estimators=args.num_trees,
            n_jobs=args.num_jobs)    
    func_name = get_function_name(method)    
    sys.stderr.write("Training %s for %s crossfold validation\n" %
            (func_name,n_folds))
    ### Need to perform Leave-one-fold-out CV to compare.
    meth_type = args.type
    if meth_type == "5mC_5hmC":
        trainData = np.concatenate((X_data,Z_data),axis=1)
    elif meth_type == "5mC":
        trainData = copy.deepcopy(X_data)
    elif meth_type == "5hmC":
        trainData = copy.deepcopy(Z_data)

    Y_pred_prob = [list() for x in range(len(Y_data))]
    Y_pred = [0]*len(Y_data)
    ### Randomly split indexes
    if args.strat_kfold:
        skf = sklearn.model_selection.StratifiedKFold(n_splits=n_folds,shuffle=True)
        skf.get_n_splits(trainData,Y_data)
        kf_split = skf.split(trainData,Y_data)
    else:
        kf = sklearn.model_selection.KFold(n_splits=n_folds,shuffle=True)
        kf.get_n_splits(trainData)
        kf_split = kf.split(trainData)
    fold = 0
    for train_idx,test_idx in kf_split:    
        fold+=1
        sys.stderr.write("\tNow on fold %s.\n" % fold)
        X_train = trainData[train_idx]
        Y_train = Y_data[train_idx]
        Z_train = Z_data[train_idx]
        I_train = I_data[train_idx]
        sys.stderr.write("\t\tnumber of training examples: %s\n" %
                (len(X_train)))
        if args.equal_class:
            sys.stderr.write("\t\tsampling to equal classes\n")
            X_train,Y_train,Z_train,I_train = \
                create_random_equivalent_training_classes(X_train,Y_train,Z_train,I_train)
            sys.stderr.write("\t\tnow there are %s training examples.\n" %
                    (len(X_train)))

        X_test = trainData[test_idx]
        Y_test = trainData[test_idx]
        method.fit(X_train,Y_train)
        Y_pred_prob_el= method.predict_proba(X_test)
        Y_pred_el = method.predict(X_test)
        for j,y_pred in enumerate(Y_pred_el):
            Y_pred_prob[test_idx[j]] = Y_pred_prob_el[j]
            Y_pred[test_idx[j]] = y_pred            
        if args.featureImportance:
            sys.stderr.write("\t\tprinting features importances\n")
            featImp = "\t".join(method.feature_importances_.astype(str))
            featImp_H = open(".".join([sample_name,func_name,meth_type,str(fold),"featureImportance"]),'w')
            featImp_H.write(featImp + "\n")
            featImp_H.close()
    #ofh = open(".".join([sample_name,func_name,meth_type,str(args.num_trees),"pred"]),'w')  
    ofh = open(".".join([sample_name,func_name,meth_type,"pred"]),'w')  
    for j,p in enumerate(Y_pred_prob):
        l = [str(x) for x in p]
        ofh.write("\t".join(l)+"\n")
    ofh.close()

    return

        

def loso_evaluation(samples,args): 
    for i,sample in enumerate(samples):
        sys.stderr.write("Preparing evaluation for %s\n" % (sample.sample_name))
        if args.obs_sample:
            train_samples = [s for j,s in enumerate(samples) if i!=j]
        else:
            train_samples = list()
            for s in samples:
                sys.stderr.write("\tchecking sample names: %s,%s\n" % (s.sample_name,sample.sample_name))
                include = True
                for idsn in s.individual_sample_names:
                    if idsn in sample.individual_sample_names:
                        include = False
                if include:
                    sys.stderr.write("\tIncluding training sample: %s\n" % (s.sample_name))
                    train_samples.append(s)        
        ## 5mC
        X_train = list()
        for ts in train_samples:
            X_train+=ts.X            
        X_train = np.asarray(X_train)
        X_test = np.asarray(sample.X)        
        ## 5hmC
        Z_train = list()
        for ts in train_samples:
            Z_train+=ts.Z            
        Z_train = np.asarray(Z_train)
        Z_test = np.asarray(sample.Z)        
        ## Expression        
        Y_train = list()
        for ts in train_samples:
            Y_train+=ts.Y
        Y_train = np.asarray(Y_train)
        #JOHN ADD: Add something to test that Y_train is >0 otherwise you get an odd error later
        Y_test = np.asarray(sample.Y)
        ## Gene ID        
        I_train = list()
        for ts in train_samples:
            I_train+=ts.I
        I_train = np.asarray(I_train)
        I_test = np.asarray(sample.I)
        sys.stderr.write("Sample: %s, Training Labels:\n" % (sample.sample_name))
        print_label_counter(Y_train)        
        sys.stderr.write("Sample: %s, Testing Labels:\n" % (sample.sample_name))
        print_label_counter(Y_test)
        if args.equal_class:
            X_train,Y_train,Z_train,I_train = create_random_equivalent_training_classes(X_train,Y_train,Z_train,I_train)
        sys.stderr.write("Sample: %s, Training Labels:\n" % (sample.sample_name))
        print_label_counter(Y_train)        
        expression_prediction_loso(X_train,Y_train,I_train,Z_train,X_test,Y_test,I_test,Z_test,sample.sample_name,args) 

def create_random_equivalent_training_classes(X,Y,Z,I):
    counter = collections.Counter(Y)
    min_count = min(counter.values())
#    print min_count
    pos_lab_idx = np.where(Y==1)[0]
#    print pos_lab_idx
    neg_lab_idx = np.where(Y==-1)[0]
#    print neg_lab_idx
    pos_lab_randsub_idx = np.random.choice(pos_lab_idx,min_count,replace=False)
    neg_lab_randsub_idx = np.random.choice(neg_lab_idx,min_count,replace=False)
    rand_sub_idx = np.sort(np.concatenate((pos_lab_randsub_idx,neg_lab_randsub_idx),axis=0))
    return (X[rand_sub_idx],Y[rand_sub_idx],Z[rand_sub_idx],I[rand_sub_idx])

def expression_prediction_loso(M_train,Y_train,I_train,Z_train,M_test,Y_test,I_test,Z_test,sample_name,args):               
    method = sklearn.ensemble.RandomForestClassifier(n_estimators=args.num_trees,n_jobs=args.num_jobs)    
    func_name = get_function_name(method)    
#    outname = ".".join([sample_name,func_name,"png"])
    sys.stderr.write("Training %s for %s\n" % (func_name,sample_name))
    ### Need to perform Leave-one-fold-out CV to compare.
    meth_type = args.type
    if meth_type == "5mC_5hmC":
        X_train = np.concatenate((M_train,Z_train),axis=1)
        X_test = np.concatenate((M_test,Z_test),axis=1)
    elif meth_type == "5mC":
        X_train = copy.deepcopy(M_train)
        X_test = copy.deepcopy(M_test)
    elif meth_type == "5hmC":
        X_train = copy.deepcopy(Z_train)
        X_test = copy.deepcopy(Z_test)
    else:
        sys.stderr.write("ERROR, Invalid type: %s. Check help for correct values.\n" %
                (meth_type))
        return
    
    n_folds = args.k_folds
    Y_pred_prob = [list() for x in range(len(Y_test))]
    Y_pred = [0]*len(Y_test)
    ### Randomly split indexes
    if args.strat_kfold:
        #kf = sklearn.cross_validation.StratifiedKFold(Y_test,n_folds=n_folds,shuffle=True)
        kf = sklearn.cross_validation.StratifiedKFold(Y_test,n_splits=n_folds,shuffle=True)
    else:
        #kf = sklearn.cross_validation.KFold(len(Y_test),n_folds=n_folds,shuffle=True)
        kf = sklearn.model_selection.KFold(len(Y_test),n_splits=n_folds,shuffle=True)
    for train_idx,test_idx in kf:    
        I_kf_test = I_test[test_idx]
        I_kf_test_set = set(I_kf_test.tolist())
        X_train_loo = list()
        Y_train_loo = list()
        ## Need to exclude any examples of exactly the gene        
        for j,i_train in enumerate(I_train):
            if args.obs_gene:
                X_train_loo.append(X_train[j])
                Y_train_loo.append(Y_train[j])
            else:
                if not i_train in I_kf_test_set:
                    X_train_loo.append(X_train[j])
                    Y_train_loo.append(Y_train[j])
        X_train_loo = np.array(X_train_loo)
        Y_train_loo = np.array(Y_train_loo)
        method.fit(X_train_loo,Y_train_loo)
        X_test_loo = X_test[test_idx]
        Y_pred_prob_el= method.predict_proba(X_test_loo)
        Y_pred_el = method.predict(X_test_loo)
        for j,y_pred in enumerate(Y_pred_el):
            Y_pred_prob[test_idx[j]] = Y_pred_prob_el[j]
            Y_pred[test_idx[j]] = y_pred            
    ofh = open(".".join([sample_name,func_name,meth_type,"pred"]),'w')  
    for j,p in enumerate(Y_pred_prob):
        l = [str(x) for x in p]
        ofh.write("\t".join(l)+"\n")
    ofh.close()
    return

def load_samples(interp_dir):
    samples = list()
    for f in os.listdir(interp_dir):
        #if f.endswith(".tss.meth.dat"):
        match = re.search("^(\S+).meth.dat$", f)
        #sys.stderr.write("Processing file %s\n" % (f))
        if match:
            sample_name = f.split(".")[0]    
            base = match.group(1)
            sys.stderr.write("Processing sample %s, base: %s\n" % 
                    (sample_name, base))    
            dat_path = os.path.join(interp_dir,f)
            label_path = os.path.join(interp_dir,sample_name+".label")
            #hmC_path = os.path.join(interp_dir,sample_name+".tss.hmC.dat")
            hmC_path = os.path.join(interp_dir,base+".hmC.dat")
            Y,L,E,I = load_multilabel_expr(label_path)
            X = load_vals(dat_path)
            Z = load_vals(hmC_path)
            sample = Sample(sample_name,X,Y,L,E,Z,I)
            samples.append(sample)    
            print_label_counter(sample.Y)
    return samples
   
    
def print_label_counter(Y):
    counter = collections.Counter(Y)
    for k in sorted(counter.keys()):
        sys.stderr.write("Class: %s, #: %d\n" % (str(k),counter[k]))
    sys.stderr.write("\n")
    return

def load_vals(filename):
    sys.stderr.write("Loading %s\n" % filename) 
    fh = open(filename,'r')
    C = list()
    for line in fh:
        if line.startswith("#"):
            continue
        ll = [float(x) for x in line.strip().split()]
        C.append(ll)
    return C
        
def load_multilabel_expr(filename):
    fh = open(filename,'r')
    Ex = list()
    L = list()
    E = list()
    Y = list()
    I = list()
    for line in fh:
        if line.startswith("#"):
            ll = line.lstrip("#").rstrip().split()
            ii = ll.index("GENE_ID")            
            ei = ll.index("EXPR")
            xi = ll.index("NUM_EXONS")
            pli = ll.index("POS_LOW")            
            phi = ll.index("POS_HIGH")
            continue
        ll = line.strip().split()
        id_name = ll[ii]
        expr = float(ll[ei])
        E.append(int(ll[xi]))
        L.append(math.log(int(ll[phi])-int(ll[pli]),10))
        Ex.append(expr)
        I.append(id_name)
        if expr > 0:
            Y.append(1)
        else:
            Y.append(-1)
    return Y,L,E,I 

def get_function_name(func):
    return str(func).split('(')[0]

def setup_parser_arguments():
    global parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,description='''
For classification of gene expression based on interpolated methylation values over multiple samples''')
    ### Required positional arguments
    parser.add_argument('interp_dir',help="Directory of <sample>.label and <sample>.wsg.meth.dat files")
    parser.add_argument('-j','--num_jobs',default=1,type=int,
            help="Number of jobs for multithreaded methods. Use -1 for all \
            cores. (default: 1)")
    parser.add_argument('-k','--k_folds',default=10,type=int,help="k for \
            K-fold cross validation with default evalution framework \
            (default: 10)")
    parser.add_argument('-t','--type',choices=['5mC_5hmC', '5mC', '5hmC'],default="5mC_5hmC",
            help="Type of analysis.  (default: 5mC_5hmC)")
    RF_group = parser.add_argument_group(title='RandomForest Options')
    RF_group.add_argument('--num_trees',default=5001,type=int,
            help="Number of trees for the Random Forest. (default: 5001)")
    RF_group.add_argument('--featureImportance',dest='featureImportance',
            action='store_true',
            help="Print file of featureImportances for each fold.")
    RF_group.set_defaults(featureImportance=False)
    advanced_group = parser.add_argument_group(title='Advanced Training Options',
            description='Warning, improperly changing these can lead to severe overfitting!')
    advanced_group.add_argument('--stratified_kFold',dest='strat_kfold',
            action='store_true',
            help="Use stratified kFold cross validation.")
    advanced_group.set_defaults(strat_kfold=False)
    advanced_group.add_argument('--no-equal_class',dest='equal_class',
            action='store_false',
            help="Do not randomly subsample the dominant class to create equal\
            classes.")
    advanced_group.set_defaults(equal_class=True)
    advanced_group.add_argument('--obs_gene',dest='obs_gene',
            action='store_true',
            help="Include genes from training set if observed in \
            testing set. Ignored for single sample case.")
    advanced_group.set_defaults(obs_gene=False)
    advanced_group.add_argument('--obs_sample',dest='obs_sample',
            action='store_true',
            help="Include sample from testing set in training set. Ignored for\
            single sample case.")
    advanced_group.set_defaults(obs_sample=False)
    
    
    return parser

## if __name__ == "__main__":
##     import sys, os, argparse, linecache, operator, math, pickle, collections, random, copy
##     import re
##     import numpy as np
##     #import mlpy
##     import multiprocessing
##     import scipy.stats
##     import sklearn.cross_validation
##     import sklearn.ensemble
##     import sklearn.metrics
##     import sklearn.neighbors    
##     import sklearn.externals
##     
## #    from matplotlib import pyplot as plt
##     main()    
    
    
    
