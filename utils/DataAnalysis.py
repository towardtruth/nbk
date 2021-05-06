"""
    Data Analysis
        - PrecisionRecallCurve

    Author: Kyoung Tak Cho
    Created: Thu Jun 27 10:44:16 CDT 2019
    Updated: Fri Apr 17 10:48:15 CDT 2020
"""

from sklearn import svm, datasets
from sklearn.model_selection import train_test_split
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from inspect import signature

from utils.FileManager import FileType, CSVLoader
from utils.Common import ListTool


class BoxPlot(object):
    """
        Draw Box Plot for data analysis and visualization
    """
    def __init__(self, src=None):
        self.src = None
        self.df = None
        self.df_filtered = None
        self.showfliers = True
        self.is_show_legend = True
        self.font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

        if src is None:
            raise ValueError("src is None. Try again.")
        else:
            self.src = src
            self.load(src=self.src)
            self.filter()   # initialization

    def load(self, src):
        """
        Load data into DataFrame type
          structure: exp_num[0], k value[1], seq_type[2], label_type[3], label_name[4], tissue_num[5],
                     accuracy[15], precision[16], recall[17], f-measure[18]
        """
        if not os.path.exists(src):
            raise ValueError("src({}) does not exist.".format(src))

        self.df = pd.read_csv(src)

        # rename columns to all lower case and replace '-' to '_'
        self._cols_rename()

    def show_legend(self, is_show_legend=True):
        self.is_show_legend = is_show_legend

    def _cols_rename(self):
        # 1. convert lower case
        self.df = self.df.rename(str.lower, axis='columns')
        new_columns = list()

        for col in self.df.columns:
            #print(col)
            new_name = col.replace('-', '_')
            new_name = new_name.replace('k size', 'k')
            new_name = new_name.replace(' ', '_')
            new_columns.append(new_name)

        # expression level (label_name): N -> T

        # expression type (label_type): gene -> mRNA abundance, protein -> protein abundance, combined abundance: RA, PA, CA

        # sequency type: k-mer amino acid, k-mer nucleotide: kA, kN


        # rename columns
        self.df.columns = new_columns

        #for col in self.df.columns:
        #    print(col)

    def set(self, showfliers=None):
        if showfliers is not None:
            self.showfliers = showfliers

    def show(self, features=None, x=None, y=None, hue=None, save_only=False, figure_file=None,
             palette=None, showfliers=None):
        """
        Show Box Plot with filtering
          filtering parameters: k vale, label type, and tissue ID
            k: k value
            seq_type: sequence type {protein|promoter}
            label_type: labeling type {gene expression|protein abundance}
            percentile: percentile models {Top 5:30:5, Bottom 5:30:5, combination of Top and Bottom}
              detail percentile 1: {T|B|C} (Top, Bottom, Combination)
              detail percentile 2: {5|10|15|20|25|30} or {5:5..5:30..30:5..30:30 combinations}
            tissue ID
          None is the default value for all parameters which means considering all data (no filtering).
        """
        #print('test: show()')
        #print(self.df.head())

        if features is not None and len(features) == 3:
            x = features[0]
            y = features[1]
            hue = features[2]
        if showfliers is not None:
            self.set(showfliers=showfliers)

        #plt.figure(figsize=(30, 15))
        plt.figure(figsize=(30, 18))
        #plt.figure(figsize=(30, 20))
        plt.rcParams.update({'font.size': 20})
        #plt.ylim(0.1, 1)
        #plt.ylim(0, 1)
        #plt.ylim(0.45, 1)
        plt.ylim(0.4, 0.8)

        #print(self.df)
        #self.df.plot.box()
        if palette is None:
            sns.set(style="ticks", palette="pastel", font_scale=2)
        else:
            sns.set(style="ticks", palette=palette, font_scale=2)
        #sns.boxplot(x=x, y=y, data=self.df)
        #ax = sns.boxplot(x=x, y=y, hue=hue, data=self.df)
        ax = sns.boxplot(x=x, y=y, hue=hue, data=self.df_filtered, showfliers=self.showfliers)
        #if hue == 'label_name':
        #    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper right', title=hue, ncol=2)
        #else:
        #    #ax.legend(loc='upper right', title=hue, ncol=1)
        #    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper right', title=hue, ncol=1)
        #ax.legend(bbox_to_anchor=(1.1, 1), loc='upper right', title=hue, ncol=1)
        #ax.legend(bbox_to_anchor=(1, 1), loc='upper right', title=hue, ncol=1)
        ax.legend(bbox_to_anchor=(1.1, 1.15), loc='upper right', title=hue, ncol=1)
        #ax.legend(loc='lower center', title=hue, ncol=6)
        #ax.legend(loc='lower right', title=hue, ncol=1)

        #sns.despine(offset=10, trim=True)
        plt.setp(ax.get_legend().get_title(), fontsize='30')
        plt.setp(ax.get_legend().get_texts(), fontsize='30')

        if not self.is_show_legend:
            ax.get_legend().remove()

        #plt.legend(prop={'size': 20})
        plt.xlabel(x, size=40)
        plt.ylabel(y, size=40)
        #plt.xticks(rotation=45, size=25)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        plt.yticks(size=40)
        plt.grid(True)

        if save_only:
            plt.savefig(figure_file)
        else:
            plt.show()
        #plt.clf()
        plt.close('all')

    def filter(self, use_df_filtered=False, query=None, col_name=None, value_list=None):
        """
        Filtering for specific data of dataframe
        input:
          col_name
          value
        output:
          self.df_filtered := set with filtered self.df
        return: None
        """
        if use_df_filtered:
            df_src = self.df_filtered
        else:
            df_src = self.df

        if query is not None:
            self.df_filtered = df_src.query(query)
        else:
            if col_name is None or value_list is None:
                self.df_filtered = df_src
            else:
                self.df_filtered = df_src[df_src[col_name] in value_list]


class PrecisionRecallCurve(object):
    """
        Precision Recall Curve
    """
    def __init__(self, sub_ax=None):
        self.data = list()
        self.y_true = list()
        self.y_score = list()
        self.p_value1 = list()
        self.p_value2 = list()
        self._p_max = 0
        self._p_min = 0
        self.average_precision = 0
        self.file_name_csv = None
        self.sub_ax = sub_ax

    def load(self, file_name):
        """
            load csv file (prediction details file)
        :param file_name:
        :return:
            set self.data (list) with csv file
        """
        # check is csv file
        if FileType.is_csv(file_name):
            self.file_name_csv = file_name
            self.data = CSVLoader.csv2list(file_name)
        else:
            raise ValueError('No csv file detected.')
        #return self.data

    def pr_curve(self):
        self.get_column_values()
        self.p_max()
        self.p_min()
        self.comp_y_score()
        self.avg_precision_score()
        self.plot_pr_curve(is_save_fig=True, is_show_plot=False)

    def p_max(self):
        max1 = max(self.p_value1)
        max2 = max(self.p_value2)
        self._p_max = max(max1, max2)
        return self._p_max

    def p_min(self):
        min1 = min(self.p_value1)
        min2 = min(self.p_value2)
        self._p_min = min(min1, min2)
        return self._p_min

    def get_column_values(self):
        self.y_true = list()
        self.p_value1 = list()
        self.p_value2 = list()
        skip_one_row = True

        for row in self.data:
            if skip_one_row:
                skip_one_row = False
            else:
                p_value1 = float(row[4])
                p_value2 = float(row[5])

                self.y_true.append(int(row[2]))
                self.p_value1.append(p_value1)
                self.p_value2.append(p_value2)

    def comp_y_score(self):
        self.y_score = list()
        skip_one_row = True

        for row in self.data:
            if skip_one_row:
                skip_one_row = False
            else:
                p_value1 = float(row[4])
                p_value2 = float(row[5])
                ls = p_value1 / p_value2
                self.y_score.append((ls - self._p_min) / (self._p_max - self._p_min))

        return self.y_score

    def avg_precision_score(self):
        self.average_precision = average_precision_score(self.y_true, self.y_score)
        print('Average precision-recall score: {0:0.2f}'.format(self.average_precision))
        return self.average_precision

    def plot_pr_curve(self, is_save_fig=True, is_show_plot=False):
        #
        # Plot the Precision-Recall curve
        #
        precision, recall, _ = precision_recall_curve(self.y_true, self.y_score)

        # In matplotlib < 1.5, plt.fill_between does not have a 'step' argument
        step_kwargs = ({'step': 'post'}
                       if 'step' in signature(plt.fill_between).parameters
                       else {})
        #plt.step(recall, precision, color='gray', alpha=1, where='post', linewidth=0.5)
        plt.step(recall, precision, alpha=1, where='post', linewidth=0.5, label='ts# {}'.format(self.sub_ax))
        #plt.fill_between(recall, precision, alpha=0.2, color='b', **step_kwargs)

        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        plt.title('2-class Precision-Recall curve: AP={0:0.2f}'.format(self.average_precision))

        #if self.sub_ax == 23:
        #leg = plt.legend(loc='upper right', ncol=2, bbox_to_anchor=(0.1, 0.1))
        #leg = plt.legend(loc='upper right', ncol=2)
        leg = plt.legend(loc='lower left', ncol=4, prop={'size': 7})
        #leg_texts = leg.get_texts()
        #plt.step(leg_texts, fontsize='x-small')

        if is_save_fig:
            plt.savefig(self.file_name_csv[:-3] + 'pdf')
        if is_show_plot:
            plt.show()

        #plt.clf()



class DistributionOfGenes(object):
    """
        Draw reverse cumulative distribution using histograms (matplotlib)
    """
    def __init__(self):
        self.data = list()


    def load(self, file_name):
        """
            load csv file (distribution daat file from DB)
        :param file_name:
        :return:
            set self.data (list) with csv file
        """
        # check is csv file
        if FileType.is_csv(file_name):
            self.file_name_csv = file_name
            self.data = CSVLoader.csv2list(file_name, skip_header=True, numeric=True)
        else:
            raise ValueError('No csv file detected.')

    def plot_gene_dist(self):
        flg, ax = plt.subplots(figsize=(20, 4))

        #data = self.data[19800:20000]
        data = self.data
        data = ListTool.twoD2oneD(data)

        # Set bins
        l_bins = list()
        l_bins.append(0)
        for i in range(0, 4000):
            l_bins.append(pow(10, i/1000))


        # plot the reverse cumulative histogram
        #n, bins, patches = ax.hist(data, bins=[0,0.001, 0.002, 0.1, 0.2, 0.3, 0.5, 0.8, 1,2,3,4,5,1000,2000], density=False, histtype='step',
        #n, bins, patches = ax.hist(data, bins=[0, 1, 2,3,4,5,6,7,8,9, 10,11,12,13,14,15,16,17,18,19, 20,30,40,50,60,70,80,90, 100, 1000, 10000], density=False, histtype='step',
        #n, bins, patches = ax.hist(data, bins='auto', density=False, histtype='step',
        n, bins, patches = ax.hist(data, bins=l_bins, density=False, histtype='step',
                                   cumulative=-1, label='Tissue #')
        print(len(n), len(bins), len(patches))
        #print(bins)

        # tidy up the figure
        ax.legend(loc='right')
        ax.set_title('Genes distribution by tissues')
        ax.set_xlabel('Cutoff value')
        ax.set_ylabel('Number of genes')
        ax.set_xscale('log')

        #plt.xscale('log')
        plt.show()


