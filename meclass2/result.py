# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 20:14:22 2018

@author: Chris
"""

class Result:
    def __init__(self,sample_name,legend,labels,predictions,pred_scores,expression,gene_IDs,gene_names,
            sample_names_in_label_file):
        self.sample_name = sample_name
        self.sample_names_in_label_file = sample_names_in_label_file
        self.legend = legend
        self.labels = labels
        self.predictions = predictions
        self.pred_scores = pred_scores
        self.expression = expression
        self.gene_IDs = gene_IDs
        self.gene_names = gene_names

        self.roc_fpr = ()
        self.roc_tpr = ()
        self.roc_thresholds = ()
        self.roc_auc = ()

        self.acc_rejrate_accuracy = ()
        self.acc_rejrate_rejectRate = ()
        self.acc_rejrate_numGenes_P = ()
        self.acc_rejrate_numGenes_N = ()
        self.max_threshold_given_accuracy = 0.

#    def set_roc(self, fpr, tpr, thresholds):
#        self.roc_fpr = fpr
#        self.roc_tpr = tpr
#        self.roc_thresholds = thresholds

    def pull_roc(self):
        result_txt = "## " + self.legend + "," + self.sample_name + "\n"
        result_txt = result_txt + "#fpr,tpr\n"
        for i, fpr in enumerate(self.roc_fpr):
            tpr = self.roc_tpr[i]
            result_txt = result_txt + "\t".join((str(fpr),str(tpr))) + "\n"
        #print(result_txt)
        return result_txt

    def pull_acc_rejectrate(self):
        result_txt = "## " + self.legend + "," + self.sample_name + "\n"
        result_txt = result_txt + "#1-rejectRate,accuracy\n"
        for i, accuracy in enumerate(self.acc_rejrate_accuracy):
            rejectRate = self.acc_rejrate_rejectRate[i]
            result_txt = result_txt + "\t".join((str(rejectRate),str(accuracy))) + "\n"
        #print(result_txt)
        return result_txt

    def pull_high_accuracy_genes(self):
        result_list = []
        for i,pred_score in enumerate(self.pred_scores):
            if pred_score >= 1-self.max_threshold_given_accuracy or pred_score <= self.max_threshold_given_accuracy:
                result_list.append("\t".join((self.sample_name,self.sample_names_in_label_file[i],self.gene_IDs[i],self.gene_names[i],str(self.expression[i]),str(self.labels[i]),str(self.predictions[i]),str(pred_score))))

        return result_list

##  
##      def find_min_pred_score_for_given_acc(self, min_acc):
##          min_pred_score = 1.
##          for i, accuracy in enumerate(self.acc_rejrate_accuracy):
##              if accuracy > min_acc and self.acc_rejrate_rejectRate[i] < min_pred_score:
##                  min_rejrate = self.acc_rejrate_rejectRate[i]
##          return min_pred_score

