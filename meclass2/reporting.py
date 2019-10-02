# -*- coding: utf-8 -*-
"""
Created on June 8th 2017

@author: cschlosberg & jedwards
"""

from meclass2.sample_reporting import Sample
from meclass2.result import Result

import sys, os, argparse, math, copy
from collections import defaultdict
import collections
import re
import numpy as np
import sklearn.metrics
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib as mpl

def set_matplotlib_params(args, version):
        # See https://matplotlib.org/users/customizing.html
        lineWidth = float(args.lineWidth)
        params = {}
        if version >= 2:
            params = {'legend.fontsize': 25,
                    'legend.frameon': False,
                    'legend.borderaxespad': 0.,
                    'legend.handlelength': 1,
                    'legend.handletextpad': 0.3,
                    'figure.figsize': (10, 10),
                    'lines.linewidth' : lineWidth,
                    'axes.linewidth' : 2,
                    'axes.labelsize': 30,
                    'axes.titlesize':'xx-large',
                    'axes.spines.top' : False,
                    'axes.spines.right' : False,
                    'xtick.direction' : 'out',
                    'ytick.direction' : 'out',
                    'xtick.top' : False,
                    'ytick.right' : False,
                    'xtick.major.size' : 15,
                    'ytick.major.size' : 15,
                    'xtick.major.width' : 2,
                    'ytick.major.width' : 2,
                    'xtick.labelsize': 25,
                    'ytick.labelsize': 25}
        else:
            warn1 = ("WARNING: Using matplot lib version < 2. Your are using " +
                "version ")
            warn2 = ("Everything should be okay, but no guarantees.\n")
            sys.stderr.write(warn1 + "%s.  " % mpl.__version__ + warn2)
            params = {'legend.fontsize': 25,
                    'legend.frameon': False,
                    'legend.borderaxespad': 0.,
                    'legend.handlelength': 1,
                    'legend.handletextpad': 0.3,
                    'figure.figsize': (10, 10),
                    'lines.linewidth' : lineWidth,
                    'axes.linewidth' : 2,
                    'axes.labelsize': 30,
                    'axes.titlesize':'xx-large',
                    'axes.spines.top' : False,
                    'axes.spines.right' : False,
                    'xtick.direction' : 'in',
                    'ytick.direction' : 'in',
                    #'xtick.top' : False,
                    #'ytick.right' : False,
                    'xtick.major.size' : 15,
                    'ytick.major.size' : 15,
                    'xtick.major.width' : 2,
                    'ytick.major.width' : 2,
                    'xtick.labelsize': 25,
                    'ytick.labelsize': 25}
        pylab.rcParams.update(params)
        return 0

def print_high_acc_genes(results, args):
    outFile = args.outFileBase + ".geneList.txt"
    sys.stderr.write("writing high accuracy (>" + str(args.min_accuracy) + ") to " + outFile + "\n")
    header="#PREDFILE\tSAMPLE\tGENE_ID\tGENE_NAME\tEXPR\tLABEL\tPRED_LABEL\tRF_SCORE\n"
    fh = open(outFile, 'w')
    fh.write(header)
    for result in results:
        for line in result.pull_high_accuracy_genes():
            fh.write(line + "\n")
    fh.close()
    return 0


def print_acc_rejectrate (results, type, args):
    outFile = args.outFileBase + "." + type + ".acc_rejectrate.txt"
    sys.stderr.write("writing Accuracy and RejectRate data to: %s\n" % outFile)
    fh = open(outFile, 'w')
    for result in results:
        fh.write(result.pull_acc_rejectrate())
    fh.close()
    return 0

def plot_acc_rejectrate (results, outFile, args,version):
    sys.stderr.write("plotting Accuracy and RejectRate data to: %s\n" % outFile)
    fig = plt.figure()
    for result in results:
        plt.plot(result.acc_rejrate_rejectRate, result.acc_rejrate_accuracy, label=result.legend)
    plt.xlabel("1 - Reject Rate")
    plt.ylabel("Accuracy")
    axes=plt.gca()
    axes.set_xlim([0.,1.])
    axes.set_ylim([0.,1.])
    axes.set_xticks(np.arange(0., 1.1, 0.1))
    axes.set_yticks(np.arange(0., 1.1, 0.1))
    if version < 2:
        axes.yaxis.set_ticks_position('left')
        axes.xaxis.set_ticks_position('bottom')
    legendLoc = args.legendLoc
    if legendLoc == "right":
        lgd=plt.legend(bbox_to_anchor=(1.1, 1),loc=2)
    else:
        lgd=plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(outFile, bbox_extra_artists=(lgd,), bbox_inches='tight')

def print_roc (results, args):
    outFile = args.outFileBase + ".roc.txt"
    sys.stderr.write("writing ROC data to: %s\n" % outFile)
    roc_fh = open(outFile, 'w')
    for result in results:
        roc_fh.write(result.pull_roc())
    roc_fh.close()
    return 0

def plot_roc (results, outFile, args, version):
    sys.stderr.write("plotting ROC data to: %s\n" % outFile)
    fig = plt.figure()
    outFile_auc = args.outFileBase + ".roc.auc.txt"
    sys.stderr.write("writing ROC AUC data to: %s\n" % outFile_auc)
    roc_auc_fh = open(outFile_auc, 'w')
    roc_auc_fh.write("Sample\tROC_AUC\n")
    for result in results:
        #plt.plot(result.roc_fpr, result.roc_tpr, label=result.legend)
        legend = result.legend
        roc_auc_fh.write(legend + "\t" + str(result.roc_auc) + "\n")
        if args.plot_auc:
            legend = legend + ", AUC=" + str(result.roc_auc)
        plt.plot(result.roc_fpr, result.roc_tpr, label=legend)
    roc_auc_fh.close()
    plt.xlabel("FPR")
    plt.ylabel("TPR")
    axes=plt.gca()
    axes.set_xlim([0.,1.])
    axes.set_ylim([0.,1.])
    axes.set_xticks(np.arange(0., 1.1, 0.1))
    axes.set_yticks(np.arange(0., 1.1, 0.1))
    axes.spines['right'].set_visible(True)
    axes.spines['top'].set_visible(True)
    if version < 2:
        axes.yaxis.set_ticks_position('left')
        axes.xaxis.set_ticks_position('bottom')
    legendLoc = args.legendLoc
    if legendLoc == "right":
        lgd=plt.legend(bbox_to_anchor=(1.1, 1),loc=2)
    else:
        lgd=plt.legend(loc="best")
    plt.savefig(outFile, bbox_extra_artists=(lgd,), bbox_inches='tight')
    return 0

def calculate_acc_rejectrate (results, type, args):
    sys.stderr.write("calculating accuracy and reject rate, class: %s\n" % type)
    steps = args.acc_rejectrate_steps
    outFile_auc = args.outFileBase + "." + type + ".acc_rr.auc.txt"
    sys.stderr.write("writing ROC AUC data to: %s\n" % outFile_auc)
    auc_fh = open(outFile_auc, 'w')
    auc_fh.write("Sample\tacc_rr_AUC\n")
    min_accuracy = args.min_accuracy
    #sys.stderr.write("\tcalculating gene list using accuracy: %s\n" % min_accuracy)
    for result in results:
        if args.verbose:
            sys.stderr.write("\ton sample %s\n" % (result.sample_name))
        accuracy_list = list()
        rejectRate_list = list()
        numGenes_P_list = list()
        numGenes_N_list = list()
        max_threshold_given_accuracy = 0.
        for threshold in np.linspace(0,0.5,steps):
            #sys.stderr.write("\t\tthresh %s\n" % (str(threshold)))
            TP_P = 0
            TP_N = 0
            numGenes_P = 0
            numGenes_N = 0
            totalGenes = 0
            totalGenes_P = 0
            totalGenes_N = 0
            for i,prediction in enumerate(result.pred_scores):
                totalGenes += 1
                if result.labels[i] < 0:
                    totalGenes_N += 1
                else:
                    totalGenes_P += 1
                if prediction >= 1-threshold:
                    numGenes_P += 1
                    if result.labels[i] == result.predictions[i]:
                        TP_P += 1
                elif prediction <= threshold:
                    numGenes_N += 1
                    if result.labels[i] == result.predictions[i]:
                        TP_N += 1
            accuracy = 0
            rejectRate = 0
            if type == 'all':
                if (numGenes_P + numGenes_N) > 0:
                    accuracy = float(TP_P + TP_N) / float(numGenes_P + numGenes_N)
                    #sys.stderr.write("\tthreshold= %s, accuracy= %s, max_threshold_given_accuracy= %s\n" % (threshold, accuracy, max_threshold_given_accuracy))
                    if accuracy > min_accuracy and threshold > max_threshold_given_accuracy:
                        max_threshold_given_accuracy = threshold
                if totalGenes > 0:
                    rejectRate = float(numGenes_P + numGenes_N) / float(totalGenes)
            elif type == "pos":
                if (numGenes_P) > 0:
                    accuracy = float(TP_P) / float(numGenes_P)
                if totalGenes_P > 0:
                    rejectRate = float(numGenes_P) / float(totalGenes_P)
            elif type == "neg":
                if (numGenes_N) > 0:
                    accuracy = float(TP_N) / float(numGenes_N)
                if totalGenes_N > 0:
                    rejectRate = float(numGenes_N) / float(totalGenes_N)
            else:
                sys.stderr.write("\tinvalid type in calculate_acc_rejectrate: %s\n" 
                        % (type))
                exit()
            accuracy_list.append(accuracy)
            rejectRate_list.append(rejectRate)
            numGenes_P_list.append(numGenes_P)
            numGenes_N_list.append(numGenes_N)
        if type == 'all':
            result.max_threshold_given_accuracy = max_threshold_given_accuracy
            sys.stderr.write("\tmax_threshold_given_accuracy= %s\n" % max_threshold_given_accuracy)
        result.acc_rejrate_accuracy = accuracy_list
        result.acc_rejrate_rejectRate = rejectRate_list
        result.acc_rejrate_numGenes_P = numGenes_P_list
        result.acc_rejrate_numGenes_N = numGenes_N_list
        acc_rr_auc = np.trapz(accuracy_list,x=rejectRate_list)
        auc_fh.write(result.legend + "\t" + str(acc_rr_auc) + "\n")
    auc_fh.close()
    return 0


def calculate_roc (results, args):
    sys.stderr.write("calculating roc\n")
    for result in results:
        if args.verbose:
            sys.stderr.write("\ton sample %s\n" % (result.sample_name))
        fpr, tpr, thresholds = sklearn.metrics.roc_curve(result.labels,
                result.pred_scores)
        auc = sklearn.metrics.roc_auc_score(result.labels, result.pred_scores)
        result.roc_fpr = fpr
        result.roc_tpr = tpr
        result.roc_thresholds = thresholds
        result.roc_auc = auc

    return 0

def calculate_labels (samples, args):
    sys.stderr.write("calculating labels\n")
    for sample in samples:
        if args.verbose:
            sys.stderr.write("\ton sample %s\n" % (sample.sample_name))
        labels = list()
        for expr in sample.expression:
            if expr > 0:
                labels.append(1)
            else:
                labels.append(-1)
        sample.set_labels(labels)
        #sys.stderr.write("\texample labels %s\n" % (str(labels[0])))
    return 0


def calculate_results(samples, expression_dict, gene_ID_dict, gene_name_dict, sample_name_dict, args):
    results = list()
    sys.stderr.write("calculating results\n")
    for sample in samples:
        if args.verbose:
            sys.stderr.write("\ton sample %s\n" % (sample.sample_name))
        predictions = list()
        for pred in sample.predictions:
            if pred >= 0.5:
                predictions.append(1)
            else:
                predictions.append(-1)
        results.append( Result(sample.sample_name,sample.legend,sample.labels,
            predictions,sample.predictions,expression_dict[sample.sample_name],gene_ID_dict[sample.sample_name],
            gene_name_dict[sample.sample_name],sample_name_dict[sample.sample_name]) )
    return results
                

def load_labels(sample_list, args):
    sys.stderr.write("loading label files\n")
    expression_data = defaultdict(list)
    gene_ID_data = defaultdict(list)
    gene_name_data = defaultdict(list)
    samples_data = defaultdict(list)
    for f in sample_list:
        if args.verbose:
            sys.stderr.write("\ton file %s\n" % (f))
        #if f.endswith(".tss.meth.dat"):
        legend = f
        sample_name = f
        match = re.search("([0-9a-zA-Z\_]+).pred$", f)
        if match:
            sample_name = f.split(".")[0]    
        if args.verbose:
            sys.stderr.write("\tSample_name: %s\n" % 
                    (sample_name))    
        labFile = ''
        if args.label_file:
            labFile = args.label_file
        else:
            labFile = sample_name + '.label'
        expression,gene_ID,gene_name,samples = load_multilabel_expr(labFile, args)
        expression_data[sample_name] = expression  
        gene_ID_data[sample_name] = gene_ID  
        gene_name_data[sample_name] = gene_name
        samples_data[sample_name] = samples
    return expression_data, gene_ID_data, gene_name_data, samples_data
   
def load_samples(sample_list, expression_data, gene_ID_data, samples_data, args):
    sys.stderr.write("loading samples\n")
    samples = list()
    for f in sample_list:
        if args.verbose:
            sys.stderr.write("\ton file %s\n" % (f))
        #if f.endswith(".tss.meth.dat"):
        legend = f
        sample_name = f
        match = re.search("([0-9a-zA-Z\_\.]+).pred$", f)
        if match:
            sample_name = f.split(".")[0]    
            legend = match.group(1)
        if args.verbose:
            sys.stderr.write("\tProcessing sample %s, legend: %s\n" % 
                    (sample_name, legend))    
        expression = expression_data[sample_name]
        gene_ID = gene_ID_data[sample_name]
        sample_names_in_label_file = samples_data[sample_name]
        predictions = load_pred(f, args)
        sample = Sample(sample_name,legend, predictions, expression,
                gene_ID, sample_names_in_label_file)
        samples.append(sample)    
    return samples
   
    
def print_label_counter(Y):
    counter = collections.Counter(Y)
    for k in sorted(counter.keys()):
        sys.stderr.write("Class: %s, #: %d\n" % (str(k),counter[k]))
    sys.stderr.write("\n")
    return

def load_pred(filename, args):
    if args.verbose:
        sys.stderr.write("\tLoading %s\n" % filename) 
    fh = open(filename,'r')
    predictions = list()
    for line in fh:
        if line.startswith("#"):
            continue
        data = line.strip().split()
        predictions.append(float(data[1]))
    return predictions
        
def load_multilabel_expr(filename, args):
    if args.verbose:
        sys.stderr.write("\tLoading %s\n" % filename) 
    fh = open(filename,'r')
    Ex = list() #expression values
    L = list()
    E = list()
    Y = list()
    I = list() #gene_IDs
    GN = list() #gene_names
    SAMP=list()
    for line in fh:
        if line.startswith("#"):
            ll = line.lstrip("#").rstrip().split()
            ii = ll.index("GENE_ID")            
            gn = ll.index("GENE_NAME")            
            ei = ll.index("EXPR")
            xi = ll.index("NUM_EXONS")
            pli = ll.index("POS_LOW")            
            phi = ll.index("POS_HIGH")
            sam = ll.index("SAMPLE")
            continue
        ll = line.strip().split()
        id_name = ll[ii]
        gene_name = ll[gn]
        expr = float(ll[ei])
        E.append(int(ll[xi]))
        L.append(math.log(int(ll[phi])-int(ll[pli]),10))
        Ex.append(expr)
        I.append(id_name)
        GN.append(gene_name)
        SAMP.append(ll[sam])
        if expr > 0:
            Y.append(1)
        else:
            Y.append(-1)
    #return Y,L,E,I 
    return  Ex,I,GN,SAMP

def get_function_name(func):
    return str(func).split('(')[0]

def setup_parser_arguments():
    global parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,description='''
Determine and plot results from ME-Class predictions.''')
    ### Required positional arguments
    parser.add_argument('pred_files',nargs='+',help="1 or more prediction files") 
    parser.add_argument('-o', '--outFileBase',dest='outFileBase',
        default='results',
        help="base tag for output files, default = \"results\"")
    parser.add_argument('--label_file',dest='label_file',default='',
            help="label file if there is a single common label file for all\
            .pred files. Otherwise (default behaviour) assume for every\
            <base>.pred file there is a corresponding <base>.label file.")
    parser.add_argument('--plot_results',dest='plot_results',
            action='store_true',
            help="Plot results.")
    parser.set_defaults(plot_results=False)
    parser.add_argument('--plot_auc',dest='plot_auc',
            action='store_true',
            help="Print ROC AUC in legend of ROC plots")
    parser.set_defaults(plot_auc=False)
    parser.add_argument('--verbose',dest='verbose',
            action='store_true',
            help="Print extra information to stderr.")
    parser.set_defaults(verbose=False)
    ### Additional Options
    additional_group = parser.add_argument_group(
            title='Additional Options')
    additional_group.add_argument('--min_accuracy',
            help="Minimum accuracy to use for gene list output. Affected by \
                acc_rejectrate_steps option below.",
            type=float,
            default=0.9)
    ### Advanced Plotting Options
    advanced_plot_group = parser.add_argument_group(
            title='Advanced Plotting Options')
    advanced_plot_group.add_argument('--acc_rejectrate_steps',
            help="Number of points to use in accuracy vs reject rate plot, 11 steps \
                    would be points every 0.1, default=101",
            type=int,
            default=101)
    advanced_plot_group.add_argument('--legendLoc',
            help="location of legend for plots, default=best",
            choices=["best","right"],
            default="best")
    advanced_plot_group.add_argument('--lineWidth',
            help="Line width for plots, default=t.5",
            type=float,
            default=2.5)
    return parser

    
    
    
