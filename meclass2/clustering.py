# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 11:30:15 2016

@author: cschlosberg & jredwards
"""
 
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

from meclass2.sample import Sample        
 
import sys, os, argparse, linecache, re
import numpy as np
import seaborn as sns
import pandas as pd
import scipy.stats
import scipy.spatial.distance
import scipy.cluster.hierarchy    
import sklearn.ensemble
import sklearn.metrics
import sklearn.cluster
import sklearn.manifold
import sklearn.neighbors    
import sklearn.externals


import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator





def process_samples(interp_dir,tag,args):    
    samples = list()
    for x in os.listdir(interp_dir):
        if x.endswith(".meth.dat"):
            sample_name = x.split(".")[0]   
            sys.stderr.write("Preparing evaluation for %s\n" %(sample_name)) 
            dat_path = os.path.join(interp_dir,x)
            label_path = os.path.join(interp_dir,sample_name+".label")

            pred_path = os.path.join(interp_dir,sample_name+".RandomForestClassifier.pred")
            if tag:
                pred_path = os.path.join(interp_dir,sample_name+".RandomForestClassifier." + tag +".pred")
            hmC_path = os.path.join(interp_dir,sample_name+".tss.hmC.dat")
            if not os.path.isfile(hmC_path):
                sys.stderr.write("\nWARNING: Can't find matching hmC dat file for \
                        %s: %s\nThis should have been generated during \
                        interpolation, if you interpolated both 5mC and \
                        5hmC.\n\n" % (x,hmC_path))
                hmC_path = dat_path
                args.cluster_data = 'meth_only'
                #callHelp(parser)
            if not os.path.isfile(pred_path):
                sys.stderr.write("\nERROR: Can't find matching pred file for %s: %s\n\n" % (x,pred_path))
                callHelp(parser)
            Y,L,E,M,gene_id_list = load_label_expr(label_path)
            X,data_info_X = load_vals(dat_path)
            C,data_info_C = load_vals(hmC_path)
            P,data_info_P = load_vals(pred_path)
            Z = calculate_ZOE(P,Y)
            sample = Sample(sample_name,X,Y,C,P,M,E,Z,gene_id_list,data_info_X) 
            samples.append(sample) 
    return samples

def convert_coord_to_bin (lb, ub, data_info):
    tssWin = data_info[0]
    tssBins = data_info[1]
    tssOffset = data_info[2]
    bin_size = (tssWin /tssBins)
    #sys.stderr.write("\tbinsize: " + str(bin_size) + "," + str(tssBins) + "\n")
    lb_shift = lb + int(tssWin/2) - int(tssOffset/2)
    ub_shift = ub + int(tssWin/2) - int(tssOffset/2)
    lb_bin = int(lb_shift / bin_size)
    ub_bin = int(ub_shift / bin_size)
    return lb_bin, ub_bin
       
def all_sample_clustering(tag,samples,lpb,upb,cp,mi_ub,mi_lb,num_clusters,args):    
    linkage_method = args.linkage_method           
    sys.stderr.write("bounds: %d, %d\n" % (mi_lb, mi_ub))
    sys.stderr.write("pred bounds: " + str(lpb) + "," + str(upb) + "\n")
    
    Y = list()
    X = list()
    Xsub = list()
    Csub = list()
    S = list()
    C = list()
    print_labels = list()
    #assume default binning scheme if problems parsing it later
    tssWin = 10000
    tssBins = 500
    tssOffset = 0
    data_info = (tssWin, tssBins, tssOffset)
    
    for s,sample in enumerate(samples):    
        idx = [i for i,z in enumerate(sample.Z) if (z==cp and (lpb<=max(sample.P[i])<=upb))]   
        sX = [x for i,x in enumerate(sample.X) if i in idx ]
        if len(sX)==0 or len(sX)==1:
            return
        if sample.data_info:
            #sys.stderr.write("\tsample.data_info:" + str(sample.data_info[0]) + "," + str(sample.data_info[1]) + "," + str(sample.data_info[2])  + "\n")
            data_info = sample.data_info
        (mi_lb_bin, mi_ub_bin) = convert_coord_to_bin( mi_lb, mi_ub, data_info)
        #sys.stderr.write("\tbins: " + str(mi_lb_bin) + "," + str(mi_ub_bin) + "\n")
        sXsub = [x[mi_lb_bin:mi_ub_bin] for i,x in enumerate(sample.X) if i in idx ]
        sY = [y for i,y in enumerate(sample.Y) if i in idx ]
    #    P = [p for i,p in enumerate(sample.P) if i in idx ]
    #    L = [math.log(l,10) for i,l in enumerate(sample.L) if i in idx ]
#        E = [e for i,e in enumerate(sample.E) if i in idx ]
        Z = [z for i,z in enumerate(sample.Z) if i in idx ]
        sC = [c for i,c in enumerate(sample.C) if i in idx ]
        #sP = [p for i,p in enumerate(sample.P) if i in idx ]
        s_labels = ["\t".join((geneID,str(sample.Y[i]),str(sample.P[i]))) for i,geneID in enumerate(sample.geneIDs) if i in idx ]
        sCsub = [c[mi_lb_bin:mi_ub_bin] for i,c in enumerate(sample.C) if i in idx ]
        sS = [s for i in range(len(sX))]
        
        Y += sY
        X += sX
        S += sS
        C += sC
        print_labels += s_labels
        
        Xsub+=sXsub
        Csub+=sCsub
    
    F = list()

    if args.cluster_data == 'meth_only':
        for i,x in enumerate(Xsub):
            F.append(x)
    elif args.cluster_data == 'hmC_only':
        for i,x in enumerate(Csub):
            F.append(x)
    else:
        for i,x in enumerate(Xsub):
            F.append(x + Csub[i])
    F = sklearn.preprocessing.normalize(F)
    
    if cp:
        sys.stderr.write("Number of genes correctly predicted with probability of classification [%f,%f]: %d in %d samples\n"%(lpb,upb,len(X),len(samples)))
    else:
        sys.stderr.write("Number of genes incorrectly predicted with probability of classification [%f,%f]: %d in %d samples\n"%(lpb,upb,len(X),len(samples)))
    of_base = "%s.lb_%f.ub_%f.cp_%d" % (tag,lpb,upb,cp)

    
### Methylation plot
    if args.cluster_data != 'hmC_only':
        title = "Differential Methylation"
        ofn_meth = of_base+".meth.clustermap.png"
        norm_Y = normalize_expression(Y) 
        meth_cmap = sns.diverging_palette(240,10,n=15,as_cmap=True)            
        sample_cmap = plt.get_cmap("gnuplot2")
        sample_tags = [float(1)/(s+1) for s in S]    
        cluster_cmap = plt.get_cmap("cool")
        #pred_cmap = plt.get_cmap("YlGnBu")
        pred_cmap = plt.get_cmap("RdYlGn_r")
        df = pd.DataFrame(X)
        vmin = -1
        vmax = 1
        sys.stderr.write("Clustering methylation and CpG density\n")
        linkage = scipy.cluster.hierarchy.linkage(F,method=linkage_method,metric='euclidean')
        fcluster = scipy.cluster.hierarchy.fcluster(linkage,num_clusters,criterion='maxclust')
        #cluster_tags = [float(1)/ct for ct in fcluster] 
        cluster_tags = [float(ct % 2) for ct in fcluster] 
        uniq_clusters = list(set(fcluster))
        row_colors = [sample_cmap(sample_tags),pred_cmap(norm_Y),cluster_cmap(cluster_tags)]      
        cluster_plot_helper(df,cluster_tags,row_colors,meth_cmap,linkage,ofn_meth,title,vmin,vmax,args)
    
### 5hmC plot    
    if args.cluster_data != "meth_only":
        title = "Differential 5hmC"
        ofn_hmC = of_base+".hmC.clustermap.png"
        vmin = -1
        vmax = 1
        hmC_cmap = sns.diverging_palette(240,10,n=15,as_cmap=True)  
        df = pd.DataFrame(C)
        #cluster_tags = [float(1)/ct for ct in fcluster]
        cluster_tags = [float(ct % 2) for ct in fcluster]
        #sys.stderr.write("\t\tct: " + str(len(fcluster)) + "\n")
        #for ct in fcluster:
            #sys.stderr.write("\t\tct: " + str(ct) + "," + str(ct % 2) + "\n")
        row_colors = [sample_cmap(sample_tags),pred_cmap(norm_Y),cluster_cmap(cluster_tags)]  
        cluster_plot_helper(df,cluster_tags,row_colors,hmC_cmap,linkage,ofn_hmC,title,vmin,vmax,args)


### Group all plots together 
### works but requires imageMagik, commented out to remove dependence
##      ofn = of_base+".clustermap.png"
##      cmd = "convert "+" ".join([" %s -gravity Center -crop 2600x2000+0+0 +repage " % (f) for f in [ofn_meth,ofn_hmC]])+" -append %s" %(ofn)
##      os.system(cmd)
##      os.system("convert %s -gravity Center -crop 2600x2000+0+0 +repage %s -gravity Center -crop 2600x2000+0+0 +repage -append %s" % (ofn_meth,ofn_hmC,ofn))
##  #    os.system("rm %s" % (ofn_meth))
##  #    os.system("rm %s" % (ofn_hmC)) 

    sns.reset_orig()
    
    
### Plot cluster averages for methylation/CpG density    
    print_individual_cluster_averages_meth_hmC(uniq_clusters,fcluster,X,Y,C,Z,of_base,print_labels,data_info,args)        
        
    return       
       
def cluster_plot_helper(df,cluster_tags,row_colors,val_cmap,linkage,ofn,title,vmin,vmax,args):
    sys.stderr.write("Plotting %s\n"%ofn)    
    #linkage_method = "complete"
    linkage_method = args.linkage_method
    sns.set(style="white")
    sns.clustermap(df, row_colors=row_colors,col_cluster = False,\
                           figsize=(35,25),  method=linkage_method, row_linkage=linkage,\
                           cmap=val_cmap, linewidths = 0,\
                           xticklabels=False,yticklabels=False,\
                           vmin = vmin, vmax = vmax)
    plt.title(title)
    plt.savefig(ofn)
    plt.close()   
    return

def check_purity(Y):
    direction = 'up-reg'
    if (len(Y) > 0):
        up_cnt = 0
        dn_cnt = 0
        for y in Y:
            if y>0:
                up_cnt += 1
            elif y < 0:
                dn_cnt += 1
        purity = up_cnt / (up_cnt + dn_cnt)
        if purity < 0.5:
            purity = 1-purity
            direction = 'dn-reg'
    else:
        purity = -1
    return purity,direction

def print_individual_cluster_averages_meth_hmC(uniq_clusters,fcluster,X,Y,C,Z,of_base,labels,data_info,args):    
    for cluster in uniq_clusters:
        idx = list()        
        for i,fc in enumerate(fcluster):
            if cluster == fc:
                idx.append(i)
        Xs = [x for i,x in enumerate(X) if i in idx]
        Ys = [y for i,y in enumerate(Y) if i in idx]
        Cs = [c for i,c in enumerate(C) if i in idx]
        Xs = np.asarray(Xs)
        Ys = np.asarray(Ys)
        Cs = np.asarray(Cs)
        purity,expression_direction = check_purity(Ys)
        sys.stderr.write("cluster: " + str(cluster) + 
                "; Expr. Dir: " + expression_direction +
                "; purity: " + str(purity) +
                "; num_genes: " + str(len(Ys)) + "\n")
        if len(Xs) >= args.min_genes_per_cluster and purity >= args.min_cluster_purity:
            cluster_labels = [l for i,l in enumerate(labels) if i in idx]
            average_print_helper_meth_cpg(Xs,Cs,str(cluster),of_base,cluster_labels,purity,expression_direction,data_info,args)
    return 0

def average_print_helper_meth_cpg(Xs,Cs,cluster,base,labels,purity,expression_direction,data_info,args):
    df = list()    
    sns.set(font_scale=1.8)
    sns.set_style("ticks")

    plt.figure(figsize=(8,5)) 
    for i,x in enumerate(Xs):      
        for j,y in enumerate(x):        
            if args.cluster_data == 'meth_only':
                df.append([j,i,y,r'$\Delta$mCG/CG'])
            elif args.cluster_data == '5hmC_only':
                df.append([j,i,Cs[i][j],r'$\Delta$hmCG/CG'])
            else:
                df.append([j,i,y,r'$\Delta$mCG/CG'])
                df.append([j,i,Cs[i][j],r'$\Delta$hmCG/CG'])
                df.append([j,i,y+Cs[i][j],r'$\Delta$mCG/CG+$\Delta$hmCG/CG'])
    df = pd.DataFrame(df)
    ci = args.confidence_interval
    #tsplot was replaced by lineplot in seaborn 0.9
    if "tsplot" in dir(sns):
        df.columns = [0,1,2,"signal"]
        ax = sns.tsplot(df,time=0,unit=1,value=2,condition="signal",ci=ci)
    else:
        df.columns = ["loc","dummy","meth","signal"]
        ax = sns.lineplot(data=df,x='loc',y='meth',hue="signal",ci=ci)

    max_x_val = Xs.shape[1]
    max_y_val = args.max_y_value
    (tssWin, tssBins, tssOffset) = data_info
    tss_ratio = float(tssWin) / (float(tssWin / 2.) + float(tssOffset))
    md_pt = int(float(max_x_val)/tss_ratio)
    true_md_pt = int(md_pt * float(tssWin) / float(tssBins))
    min_label =  true_md_pt - tssWin
    max_label = tssWin - true_md_pt

    plt.plot([md_pt,md_pt],[-1,1],'k-',alpha=0.5)
    plt.plot([0,max_x_val],[0,0],'k--',alpha=0.5)    
    #print(max_x_val,max_y_val,min_label,max_label,md_pt,true_md_pt)

    plt.xlim([0,max_x_val])
    #plt.xticks((0,md_pt,max_x_val),('-5kb','TSS','+5kb'))
    plt.xticks((0,md_pt,max_x_val),(min_label,'TSS',max_label))
    plt.ylim([-max_y_val,max_y_val])
    plt.yticks(np.arange(-max_y_val, max_y_val*1.05, step=0.2))
    #ax.minorticks_on()
    #x_tick_minorLocator = MultipleLocator(1000)
    #ax.xaxis.set_minor_locator(x_tick_minorLocator)

    plt.title("Cluster: %s (n=%d, %s, purity=%5.2f)" % (cluster,len(Xs),expression_direction,purity))
    #ax.legend_.remove()
    lgd=plt.legend(loc=5,bbox_to_anchor=(1.7,0.5),handlelength=1,handletextpad=0.5)
    plt.xlabel("Position relative to TSS (bp)")
    if args.cluster_data == 'meth_only':
        plt.ylabel(r'$\Delta$mCG/CG')
    elif args.cluster_data == 'hmC_only':
        plt.ylabel(r'$\Delta$hmCG/CG')
    else:
        plt.ylabel(r'$\Delta$mCG/CG'+'\n'+r'$\Delta$hmCG/CG')

    sys.stderr.write("\tPlotting Aggregate 5mC/5hmC Cluster %s\n" % (cluster))
    if args.tight_layout:
        plt.tight_layout()
    plt.savefig(base+".meth_cpg.cluster_%s"%(cluster)+".png", bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()

    outFile = (base+".meth_cpg.cluster_%s"%(cluster)+".txt")
    outFile_fh = open(outFile, 'w')
    header = "GeneID\tExpr\tPred\n"
    outFile_fh.write(header)
    gene_output_text = "\n".join(labels)
    outFile_fh.write(gene_output_text)
    outFile_fh.close()

    return 0

def count_vals(Y,val):
    sum = 0
    for y in Y:
        if y == val:
            sum+=1
    return sum


def normalize_cpg_density(C):
    cm_norm = mpl.colors.Normalize(vmin=0.,vmax=0.3)
    return cm_norm(C)

def normalize_expression(E):
    cm_norm = mpl.colors.Normalize(vmin=-1.,vmax=1.) 
    return cm_norm(E)

def calculate_accuracy(A):
    corr = A.count(1)
    return (float(corr)/len(A))*100
    
def load_methylation(index,filename):
    line = linecache.getline(filename,index)
    return [float(x) for x in line.strip().split()]


    
def calculate_ZOE(pred_prob,Y):
    Z = list()
    for i,p in enumerate(pred_prob):
        if p[0] > p[1]:
            res = -1
        else:
            res = 1
        if np.sign(Y[i]) == res:
            Z.append(1)
        else:
            Z.append(0)
    return Z
    
def make_violinplots(data,acc,cutoffs,name,ofn):
    sys.stderr.write("Processing %s\n"%ofn)    
    d = {"LB_cutoff":list(),name:list()}    
    str_acc = ["%.2f\n(%.2f%%)"%(cutoffs[i],a) for i,a in enumerate(acc)]
    plt.figure(figsize=(20,10))     
    for i,c in enumerate(cutoffs):
        if i == len(cutoffs)-1:
            continue
        for el in data[i]:
            d["LB_cutoff"].append(c)
            d[name].append(el)
    df = pd.DataFrame(d)
    sns.set(style="white",font_scale=1.5)
    vplot = sns.violinplot(x="LB_cutoff",y=name,data=df,scale="count")
    vplot.set_xticklabels(str_acc)
    vplot.set_xlabel("Probability of Prediction\n(Accuracy of Predictions)")
    plt.savefig(ofn)
    plt.close()
    
    
def load_vals(filename):
    fh = open(filename,'r')
    X = list()            
    sys.stderr.write("Loading %s\n" % filename)  
    tssWin = 10000
    tssBins = 500
    tssOffset = 0
    for line in fh:
        if line.startswith("#"):
            #tssWin:10000,tssBins:500,tssOffset:0
            match = re.search("tssWin:(\d+)",line)
            if match:
                tssWin = int(match.group(1))
            match = re.search("tssBins:(\d+)",line)
            if match:
                tssBins = int(match.group(1))
            match = re.search("tssOffset:(\d+)",line)
            if match:
                tssOffset = int(match.group(1))
            continue        
        X.append([float(x) for x in line.strip().split()])
    #sys.stderr.write( "tssBins: " + str(tssBins) + "\n" )
    data_info = (tssWin, tssBins, tssOffset)
    return X,data_info  
    
def load_label_expr(filename):
    fh = open(filename,'r')
    Y = list()
    L = list()
    E = list()
    M = list()
    gene_id_list = list()
    for line in fh:
        if line.startswith("#"):
            ll = line.lstrip("#").rstrip().split()
            ei = ll.index("EXPR")
            phi = ll.index("POS_HIGH")
            pli = ll.index("POS_LOW")            
            nei = ll.index("NUM_EXONS")
            gene_id_index = ll.index("GENE_ID")
            continue
        ll = line.strip().split()
        expr = float(ll[ei])
        length = int(ll[phi])-int(ll[pli])
        num_exons = int(ll[nei])
        gene_id = ll[gene_id_index]

        c = get_expression_class(expr)
        Y.append(expr)
        M.append(c)
        L.append(length)
        E.append(num_exons)
        gene_id_list.append(gene_id)
    return (Y,L,E,M,gene_id_list)  

def get_expression_class(expr):
    if expr < 0:
        c = -1
    else:
        c = 1    
    return c

def get_function_name(func):
    return str(func).split('(')[0]

def callHelp(parser):
    parser.print_usage()
    sys.exit()

def setup_parser_arguments():
    global parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,description='''
For plotting of gene signatures of classification of gene expression based on
interpolated methylation values over multiple samples. This script assumes that
you have interpolated both 5mC and 5hmC curves.  I.e. there need to be both
<base>.tss.meth.dat and <base>.tss.hmC.dat files corresponding to the pred file
in the interp_dir used as input.''')
    ### Required positional arguments
    parser.add_argument('interp_dir',help="Directory of <sample>.label, <sample>.wsg.meth.dat, <sample>.wsg.cpg.dat, <sample>.<method>.pred files")
    parser.add_argument('base',help="Base for output file names")
    parser.add_argument('--tag',
        help="tag to subset .pred files, e.g. if <>.RandomForestClassifier.5mC_5hmC.pred, then use tag=5mC_5hmC")
    parser.add_argument('--cluster_data', default="both", choices=["meth_only", "hmC_only", "both"],
        help="which signatures/data to cluster based on. (default: both)")
    parser.add_argument('--linkage_method', default="complete", choices=["single", "complete", "average", "weighted", "ward", "median", "centroid"],
        help="linkage method for clustering. See scipy.cluster.hierarchy.linkage online documentation for more info. (default: complete)")
    parser.add_argument('--upperBound',default=2500,type=int,
        help="upperBound of window for clustering in number of features, in bp relative to TSS (default: 2500)")
    parser.add_argument('--lowerBound',default=-500,type=int,
        help="lowerBound of window for clustering in number of features, in bp relative to TSS (default: -500)")
    parser.add_argument('--upperPredBound',default=1.0,type=float,
        help="upper prediction score bound for clustering. (default: 1.0)")
    parser.add_argument('--lowerPredBound',default=0.7,type=float,
        help="lower prediction score bound for clustering. (default: 0.7)")
    parser.add_argument('--numClusters',default=3,type=int,
        help="number of clusters. (default: 3)")
    parser.add_argument('--confidence_interval',default=95, type=int,
        help="confidence interval [0,100] for aggregate cluster plots. (default: 95)")
    parser.add_argument('--max_y_value',default=0.4, type=float,
        help="max y-value for plots. (default: 0.4)")
    parser.add_argument('--min_genes_per_cluster',default=10, type=int,
        help="min number of genes needed in a cluster to print it (default: 10)")
    parser.add_argument('--min_cluster_purity',default=0.75, type=float,
        help="min purity of the cluster to print it (default: 0.75)")
    parser.add_argument('--tight_layout', action='store_true',
        help="May help with layout of cluster plots if the plot does not fill the figure (default: False)")
    return parser
    
