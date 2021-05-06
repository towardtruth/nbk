"""
    Draw Plots
    Author: Kyoung Tak Cho
    Created: Monday, April 20, 2020 3:51:17 PM
    Updated: Monday, April 20, 2020 3:51:17 PM
"""
import settings
import os
import sys
import datetime
from utils.FileManager import FileList
from utils.DataAnalysis import BoxPlot
import argparse


plot_type_nbk = {'nbk', 'n1'}
plot_type_nbk2 = {'nbk2', 'n2'}
plot_type_nbk3 = {'nbk3', 'n3'}
plot_type_tnbk = {'tnbk', 't'}
plot_type_rda = {'rda', 'ra', 'r'}
choice_yes = {'yes', 'y', 'Yes', 'Y'}
choice_no = {'no', 'n', 'No', 'N'}

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--data_type', nargs='+',
                    choices=plot_type_nbk|plot_type_rda|plot_type_tnbk|plot_type_nbk2|plot_type_nbk3,
                    help='''Choose data type. ex) -t=nbk|tnbk|rda''')
parser.add_argument('-p', '--print_only_filters', choices=choice_yes|choice_no,
                    help='''Print only filters setting informaiton.''')

def main(argv):

    args = parser.parse_args(argv[1:])
    plot_type = None
    if len(argv) <= 1:
        parser.parse_args(['--help'])
        return

    if args.data_type:
        if set(args.data_type) & plot_type_nbk:
            plot_type = 'n'
        elif set(args.data_type) & plot_type_nbk2:
            plot_type = 'n2'
        elif set(args.data_type) & plot_type_nbk3:
            plot_type = 'n3'
        elif set(args.data_type) & plot_type_tnbk:
            plot_type = 't'
        elif set(args.data_type) & plot_type_rda:
            plot_type = 'r'
        else:
            plot_type = 'n'
    else:
        raise ValueError('No plot_type.')
        return

    save_box_plot(plot_type)


def save_box_plot(plot_type=None):
    exp_list = None
    features = None
    exp_list_nbk = ('0_3_EXP-1',
                    '0_3_EXP-2',
                    '0_3_EXP-3',
                    '0_3_EXP-4',
                    '0_3_EXP-5',
                    '0_3_EXP-6',
                    '0_3_EXP-7',
                    '0_3_EXP-8',
                   )
    exp_list_nbk_1_2 = ('0_3_EXP-1-2',)
    exp_list_nbk_1_3 = ('0_3_EXP-1-3',
                        '0_3_EXP-2-3',
                        '0_3_EXP-5-3',
                        '0_3_EXP-6-3')
    exp_list_nbk_1_4_0 = ('0_3_EXP-5_6',)
    exp_list_nbk_1_4 = ('0_4_EXP-1_2_3_4',
                        '0_4_EXP-5_6_7_8',)
    exp_list_nbk_1_4_1 = ('0_4_EXP-1_2_3_4',)
    exp_list_nbk_1_4_2 = ('0_4_EXP-5_6_7_8',)
    exp_list_nbk_1_4_3 = ('0_4_EXP-ALL',)
    exp_list_nbk_final = ('00_EXP-ALL-Final',)
    exp_list_kmd = ('00_EXP-KD-Final',)
    exp_list_ra = ('00_EXP-RA-Final',)

    exp_list_rda = ('0_RA-01',
                    '0_RA-02',
                   )

    features_1 = [('tissue_num', 'f_measure', 'label_type'),
                  ('tissue_num', 'f_measure', 'label_name'),
                  ('label_type', 'f_measure', 'label_name'),
                  ('tissue_num', 'f_measure', 'k'),
                  ('label_name', 'f_measure', 'k'),
                  ('label_name', 'f_measure', 'label_type'),
                  ('label_type', 'f_measure', 'label_type'),
                  ('k', 'f_measure', 'label_type'),
                  ]
    features_1_2 = [('tissue_num', 'f_measure', 'label_type'),
                    ('tissue_num', 'f_measure', 'label_name'),
                    ('tissue_num', 'f_measure', 'k'),
                    ('label_name', 'f_measure', 'k'),
                    ('label_name', 'f_measure', 'label_type'),
                    ]
    features_2 = [('tissue_num', 'f_measure', 'label_type'),
                    ('alphabet_size', 'f_measure', 'k'),
                    ('alphabet_size', 'f_measure', 'label_type'),
                    ('alphabet_size', 'f_measure', 'label_name'),
                    ('label_name', 'f_measure', 'alphabet_size'),
                    ('alphabet_size', 'f_measure', 'tissue_num'),
                    ('k', 'f_measure', 'alphabet_size'),
                    ('tissue_num', 'f_measure', 'alphabet_size'),
                    ('tissue_num', 'f_measure', 'label_name'),
                    ('tissue_num', 'f_measure', 'k'),
                    ('label_name', 'f_measure', 'k'),
                    ('label_name', 'f_measure', 'label_type'),
                    ]
    features_3_1 = [('tissue_num', 'f_measure', 'label_type'),
                  ('tissue_num', 'f_measure', 'label_name'),
                  ('tissue_num', 'f_measure', 'k'),
                  ('tissue_num', 'f_measure', 'seq_type'),
                  ('label_type', 'f_measure', 'tissue_num'),
                  ('label_type', 'f_measure', 'label_name'),
                  ('label_type', 'f_measure', 'k'),
                  ('label_type', 'f_measure', 'seq_type'),
                  ('label_name', 'f_measure', 'tissue_num'),
                  ('label_name', 'f_measure', 'label_type'),
                  ('label_name', 'f_measure', 'k'),
                  ('label_name', 'f_measure', 'seq_type'),
                  ('k', 'f_measure', 'tissue_num'),
                  ('k', 'f_measure', 'label_type'),
                  ('k', 'f_measure', 'label_name'),
                  ('k', 'f_measure', 'seq_type'),
                  ('seq_type', 'f_measure', 'tissue_num'),
                  ('seq_type', 'f_measure', 'label_type'),
                  ('seq_type', 'f_measure', 'label_name'),
                  ('seq_type', 'f_measure', 'k'),
                  ]
    features_3_2 = [('tissue_num', 'f_measure', 'label_type'),
                    ('tissue_num', 'f_measure', 'label_name'),
                    ('tissue_num', 'f_measure', 'algorithm'),
                    ('tissue_num', 'f_measure', 'k'),
                    ('tissue_num', 'f_measure', 'seq_type'),
                    ('label_type', 'f_measure', 'tissue_num'),
                    ('label_type', 'f_measure', 'label_name'),
                    ('label_type', 'f_measure', 'algorithm'),
                    ('label_type', 'f_measure', 'k'),
                    ('label_type', 'f_measure', 'seq_type'),
                    ('label_name', 'f_measure', 'tissue_num'),
                    ('label_name', 'f_measure', 'label_type'),
                    ('label_name', 'f_measure', 'algorithm'),
                    ('label_name', 'f_measure', 'k'),
                    ('label_name', 'f_measure', 'seq_type'),
                    ('algorithm', 'f_measure', 'tissue_num'),
                    ('algorithm', 'f_measure', 'label_type'),
                    ('algorithm', 'f_measure', 'label_name'),
                    ('algorithm', 'f_measure', 'k'),
                    ('algorithm', 'f_measure', 'seq_type'),
                    ('k', 'f_measure', 'tissue_num'),
                    ('k', 'f_measure', 'label_type'),
                    ('k', 'f_measure', 'label_name'),
                    ('k', 'f_measure', 'algorithm'),
                    ('k', 'f_measure', 'seq_type'),
                    ('seq_type', 'f_measure', 'tissue_num'),
                    ('seq_type', 'f_measure', 'label_type'),
                    ('seq_type', 'f_measure', 'label_name'),
                    ('seq_type', 'f_measure', 'algorithm'),
                    ('seq_type', 'f_measure', 'k'),
                    ]
    features_4 = [('tissue_num', 'f_measure', 'label_type'),
                  ('tissue_num', 'f_measure', 'label_name'),
                  ('tissue_num', 'f_measure', 'algorithm'),
                  ('tissue_num', 'f_measure', 'k'),
                  ('tissue_num', 'f_measure', 'seq_type'),
                  ('tissue_num', 'f_measure', 'model'),
                  ('label_type', 'f_measure', 'tissue_num'),
                  ('label_type', 'f_measure', 'label_name'),
                  ('label_type', 'f_measure', 'algorithm'),
                  ('label_type', 'f_measure', 'k'),
                  ('label_type', 'f_measure', 'seq_type'),
                  ('label_type', 'f_measure', 'model'),
                  ('label_name', 'f_measure', 'tissue_num'),
                  ('label_name', 'f_measure', 'label_type'),
                  ('label_name', 'f_measure', 'algorithm'),
                  ('label_name', 'f_measure', 'k'),
                  ('label_name', 'f_measure', 'seq_type'),
                  ('label_name', 'f_measure', 'model'),
                  ('algorithm', 'f_measure', 'tissue_num'),
                  ('algorithm', 'f_measure', 'label_type'),
                  ('algorithm', 'f_measure', 'label_name'),
                  ('algorithm', 'f_measure', 'k'),
                  ('algorithm', 'f_measure', 'seq_type'),
                  ('algorithm', 'f_measure', 'model'),
                  ('k', 'f_measure', 'tissue_num'),
                  ('k', 'f_measure', 'label_type'),
                  ('k', 'f_measure', 'label_name'),
                  ('k', 'f_measure', 'algorithm'),
                  ('k', 'f_measure', 'seq_type'),
                  ('k', 'f_measure', 'model'),
                  ('seq_type', 'f_measure', 'tissue_num'),
                  ('seq_type', 'f_measure', 'label_type'),
                  ('seq_type', 'f_measure', 'label_name'),
                  ('seq_type', 'f_measure', 'algorithm'),
                  ('seq_type', 'f_measure', 'k'),
                  ('seq_type', 'f_measure', 'model'),
                  ('model', 'f_measure', 'tissue_num'),
                  ('model', 'f_measure', 'label_type'),
                  ('model', 'f_measure', 'label_name'),
                  ('model', 'f_measure', 'algorithm'),
                  ('model', 'f_measure', 'k'),
                  ('model', 'f_measure', 'seq_type'),
                  ]
    features_final = [('tissue_num', 'f_measure', 'expression_type'),
                      ('tissue_num', 'f_measure', 'expression_level'),
                      ('tissue_num', 'f_measure', 'algorithm'),
                      ('tissue_num', 'f_measure', 'k'),
                      ('tissue_num', 'f_measure', 'sequence_type'),
                      ('tissue_num', 'f_measure', 'model'),
                      ('expression_type', 'f_measure', 'tissue_num'),
                      ('expression_type', 'f_measure', 'expression_level'),
                      ('expression_type', 'f_measure', 'algorithm'),
                      ('expression_type', 'f_measure', 'k'),
                      ('expression_type', 'f_measure', 'sequence_type'),
                      ('expression_type', 'f_measure', 'model'),
                      ('expression_level', 'f_measure', 'tissue_num'),
                      ('expression_level', 'f_measure', 'expression_type'),
                      ('expression_level', 'f_measure', 'algorithm'),
                      ('expression_level', 'f_measure', 'k'),
                      ('expression_level', 'f_measure', 'sequence_type'),
                      ('expression_level', 'f_measure', 'model'),
                      ('algorithm', 'f_measure', 'tissue_num'),
                      ('algorithm', 'f_measure', 'expression_type'),
                      ('algorithm', 'f_measure', 'expression_level'),
                      ('algorithm', 'f_measure', 'k'),
                      ('algorithm', 'f_measure', 'sequence_type'),
                      ('algorithm', 'f_measure', 'model'),
                      ('k', 'f_measure', 'tissue_num'),
                      ('k', 'f_measure', 'expression_type'),
                      ('k', 'f_measure', 'expression_level'),
                      ('k', 'f_measure', 'algorithm'),
                      ('k', 'f_measure', 'sequence_type'),
                      ('k', 'f_measure', 'model'),
                      ('sequence_type', 'f_measure', 'tissue_num'),
                      ('sequence_type', 'f_measure', 'expression_type'),
                      ('sequence_type', 'f_measure', 'expression_level'),
                      ('sequence_type', 'f_measure', 'algorithm'),
                      ('sequence_type', 'f_measure', 'k'),
                      ('sequence_type', 'f_measure', 'model'),
                      ('model', 'f_measure', 'tissue_num'),
                      ('model', 'f_measure', 'expression_type'),
                      ('model', 'f_measure', 'expression_level'),
                      ('model', 'f_measure', 'algorithm'),
                      ('model', 'f_measure', 'k'),
                      ('model', 'f_measure', 'sequence_type'),
                      ]
    features_kmd = [('tissue_num', 'f_measure', 'expression_type'),
                    ('tissue_num', 'f_measure', 'expression_level'),
                    ('tissue_num', 'f_measure', 'algorithm'),
                    ('tissue_num', 'f_measure', 'k'),
                    ('tissue_num', 'f_measure', 'km_set'),
                    ('tissue_num', 'f_measure', 'model'),
                    ('expression_type', 'f_measure', 'tissue_num'),
                    ('expression_type', 'f_measure', 'expression_level'),
                    ('expression_type', 'f_measure', 'algorithm'),
                    ('expression_type', 'f_measure', 'k'),
                    ('expression_type', 'f_measure', 'km_set'),
                    ('expression_type', 'f_measure', 'model'),
                    ('expression_level', 'f_measure', 'tissue_num'),
                    ('expression_level', 'f_measure', 'expression_type'),
                    ('expression_level', 'f_measure', 'algorithm'),
                    ('expression_level', 'f_measure', 'k'),
                    ('expression_level', 'f_measure', 'km_set'),
                    ('expression_level', 'f_measure', 'model'),
                    ('algorithm', 'f_measure', 'tissue_num'),
                    ('algorithm', 'f_measure', 'expression_type'),
                    ('algorithm', 'f_measure', 'expression_level'),
                    ('algorithm', 'f_measure', 'k'),
                    ('algorithm', 'f_measure', 'km_set'),
                    ('algorithm', 'f_measure', 'model'),
                    ('k', 'f_measure', 'tissue_num'),
                    ('k', 'f_measure', 'expression_type'),
                    ('k', 'f_measure', 'expression_level'),
                    ('k', 'f_measure', 'algorithm'),
                    ('k', 'f_measure', 'km_set'),
                    ('k', 'f_measure', 'model'),
                    ('km_set', 'f_measure', 'tissue_num'),
                    ('km_set', 'f_measure', 'expression_type'),
                    ('km_set', 'f_measure', 'expression_level'),
                    ('km_set', 'f_measure', 'algorithm'),
                    ('km_set', 'f_measure', 'k'),
                    ('km_set', 'f_measure', 'model'),
                    ('model', 'f_measure', 'tissue_num'),
                    ('model', 'f_measure', 'expression_type'),
                    ('model', 'f_measure', 'expression_level'),
                    ('model', 'f_measure', 'algorithm'),
                    ('model', 'f_measure', 'k'),
                    ('model', 'f_measure', 'km_set'),
                    ]
    features_ra =  [('tissue_num', 'f_measure', 'expression_type'),
                    ('tissue_num', 'f_measure', 'expression_level'),
                    ('tissue_num', 'f_measure', 'algorithm'),
                    ('tissue_num', 'f_measure', 'k'),
                    ('tissue_num', 'f_measure', 'grouping'),
                    ('tissue_num', 'f_measure', 'alphabet'),
                    ('expression_type', 'f_measure', 'tissue_num'),
                    ('expression_type', 'f_measure', 'expression_level'),
                    ('expression_type', 'f_measure', 'algorithm'),
                    ('expression_type', 'f_measure', 'k'),
                    ('expression_type', 'f_measure', 'grouping'),
                    ('expression_type', 'f_measure', 'alphabet'),
                    ('expression_level', 'f_measure', 'tissue_num'),
                    ('expression_level', 'f_measure', 'expression_type'),
                    ('expression_level', 'f_measure', 'algorithm'),
                    ('expression_level', 'f_measure', 'k'),
                    ('expression_level', 'f_measure', 'grouping'),
                    ('expression_level', 'f_measure', 'alphabet'),
                    ('algorithm', 'f_measure', 'tissue_num'),
                    ('algorithm', 'f_measure', 'expression_type'),
                    ('algorithm', 'f_measure', 'expression_level'),
                    ('algorithm', 'f_measure', 'k'),
                    ('algorithm', 'f_measure', 'grouping'),
                    ('algorithm', 'f_measure', 'alphabet'),
                    ('k', 'f_measure', 'tissue_num'),
                    ('k', 'f_measure', 'expression_type'),
                    ('k', 'f_measure', 'expression_level'),
                    ('k', 'f_measure', 'algorithm'),
                    ('k', 'f_measure', 'grouping'),
                    ('k', 'f_measure', 'alphabet'),
                    ('grouping', 'f_measure', 'tissue_num'),
                    ('grouping', 'f_measure', 'expression_type'),
                    ('grouping', 'f_measure', 'expression_level'),
                    ('grouping', 'f_measure', 'algorithm'),
                    ('grouping', 'f_measure', 'k'),
                    ('grouping', 'f_measure', 'alphabet'),
                    ('alphabet', 'f_measure', 'tissue_num'),
                    ('alphabet', 'f_measure', 'expression_type'),
                    ('alphabet', 'f_measure', 'expression_level'),
                    ('alphabet', 'f_measure', 'algorithm'),
                    ('alphabet', 'f_measure', 'k'),
                    ('alphabet', 'f_measure', 'grouping'),
                    ]

    #filters = ('all',
    #filters = (
    #    'all',
    filters = (
               'all',
               'seq_type in ["p"]',
               'seq_type in ["m1"]',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"]',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==4',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==5',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==7',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==3',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==4',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==7',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3 & label_type=="g"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3 & label_type=="p"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==4 & label_type=="g"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==4 & label_type=="p"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==5 & label_type=="g"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==5 & label_type=="p"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6 & label_type=="g"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6 & label_type=="p"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==7 & label_type=="g"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==7 & label_type=="p"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==3 & label_type=="g"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==3 & label_type=="p"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==4 & label_type=="g"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==4 & label_type=="p"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5 & label_type=="g"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5 & label_type=="p"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6 & label_type=="g"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6 & label_type=="p"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==7 & label_type=="g"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==7 & label_type=="p"',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==3',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==4',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==5',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==6',
               'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==7',
    )
    filters_nbk_bn = (
        'algorithm in ["NBK", "BayesNet"]',
        'algorithm in ["NBK", "BayesNet"] & seq_type in ["p"]',
        'algorithm in ["NBK", "BayesNet"] & seq_type in ["m1"]',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"]',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==4',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==5',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==7',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==3',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==4',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==7',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3 & label_type=="g"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3 & label_type=="p"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==4 & label_type=="g"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==4 & label_type=="p"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==5 & label_type=="g"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==5 & label_type=="p"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6 & label_type=="g"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6 & label_type=="p"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==7 & label_type=="g"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==7 & label_type=="p"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==3 & label_type=="g"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==3 & label_type=="p"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==4 & label_type=="g"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==4 & label_type=="p"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5 & label_type=="g"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5 & label_type=="p"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6 & label_type=="g"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6 & label_type=="p"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==7 & label_type=="g"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==7 & label_type=="p"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==3',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==4',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==5',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==6',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==7',
    )
    filters_nbk = (
        'algorithm in ["NBK"]',
        'algorithm in ["NBK"] & seq_type in ["p"]',
        'algorithm in ["NBK"] & seq_type in ["m1"]',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"]',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==4',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==5',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==7',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==3',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==4',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==7',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3 & label_type=="g"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3 & label_type=="p"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==4 & label_type=="g"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==4 & label_type=="p"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==5 & label_type=="g"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==5 & label_type=="p"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6 & label_type=="g"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6 & label_type=="p"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==7 & label_type=="g"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==7 & label_type=="p"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==3 & label_type=="g"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==3 & label_type=="p"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==4 & label_type=="g"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==4 & label_type=="p"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5 & label_type=="g"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5 & label_type=="p"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6 & label_type=="g"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6 & label_type=="p"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==7 & label_type=="g"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==7 & label_type=="p"',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==3',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==4',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==5',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==6',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==7',
    )
    filters_bn = (
        'algorithm in ["BayesNet"]',
        'algorithm in ["BayesNet"] & seq_type in ["p"]',
        'algorithm in ["BayesNet"] & seq_type in ["m1"]',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"]',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==4',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==5',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==7',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==3',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==4',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==7',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3 & label_type=="g"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3 & label_type=="p"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==4 & label_type=="g"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==4 & label_type=="p"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==5 & label_type=="g"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==5 & label_type=="p"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6 & label_type=="g"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6 & label_type=="p"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==7 & label_type=="g"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==7 & label_type=="p"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==3 & label_type=="g"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==3 & label_type=="p"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==4 & label_type=="g"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==4 & label_type=="p"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5 & label_type=="g"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5 & label_type=="p"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6 & label_type=="g"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6 & label_type=="p"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==7 & label_type=="g"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==7 & label_type=="p"',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==3',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==4',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==5',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==6',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & k==7',
    )
    filters_exps = (
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p"',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1"',
    )
    filters_exps_2 = (
        'label_name in ["N95"] & seq_type=="p"',
        'label_name in ["N90"] & seq_type=="p"',
        'label_name in ["N85"] & seq_type=="p"',
        'label_name in ["N80"] & seq_type=="p"',
        'label_name in ["N75"] & seq_type=="p"',
        'label_name in ["N70"] & seq_type=="p"',
        'label_name in ["N95"] & seq_type=="m1"',
        'label_name in ["N90"] & seq_type=="m1"',
        'label_name in ["N85"] & seq_type=="m1"',
        'label_name in ["N80"] & seq_type=="m1"',
        'label_name in ["N75"] & seq_type=="m1"',
        'label_name in ["N70"] & seq_type=="m1"',
    )
    filters_exps_3 = (
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & label_type=="p"',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & label_type=="g"',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & label_type=="p"',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & label_type=="g"',
    )

    filters_exps_4 = (
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70", "HG70_HP70", "HG60_HP60", "HG50_HP50"]',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70", "HG70_HP70", "HG60_HP60", "HG50_HP50"] & seq_type=="p"',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70", "HG70_HP70", "HG60_HP60", "HG50_HP50"] & seq_type=="m1"',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70", "HG70_HP70", "HG60_HP60", "HG50_HP50"] & seq_type=="p" & k==3',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70", "HG70_HP70", "HG60_HP60", "HG50_HP50"] & seq_type=="p" & k==4',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70", "HG70_HP70", "HG60_HP60", "HG50_HP50"] & seq_type=="p" & k==5',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70", "HG70_HP70", "HG60_HP60", "HG50_HP50"] & seq_type=="p" & k==6',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70", "HG70_HP70", "HG60_HP60", "HG50_HP50"] & seq_type=="p" & k==7',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70", "HG70_HP70", "HG60_HP60", "HG50_HP50"] & seq_type=="m1" & k==3',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70", "HG70_HP70", "HG60_HP60", "HG50_HP50"] & seq_type=="m1" & k==4',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70", "HG70_HP70", "HG60_HP60", "HG50_HP50"] & seq_type=="m1" & k==5',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70", "HG70_HP70", "HG60_HP60", "HG50_HP50"] & seq_type=="m1" & k==6',
        'label_name in ["N95", "N90", "N85", "N80", "N75", "N70", "HG70_HP70", "HG60_HP60", "HG50_HP50"] & seq_type=="m1" & k==7',
    )
    filters_exps_4_1 = (
        'label_name in ["B05"]',
        'label_name in ["B05"] & seq_type=="p"',
    )
    filters_exps_4_2 = (
        'label_name in ["N95"]',
        'label_name in ["N95"] & seq_type=="p"',
        'label_name in ["N95"] & seq_type=="m1"',
    )
    filters_exps_5 = (
        'label_name in ["N95"]',
        'label_name in ["N95"] & seq_type=="p"',
        'label_name in ["N95"] & seq_type=="m1"',
        'label_name in ["N95"] & seq_type=="p" & k==3',
        'label_name in ["N95"] & seq_type=="p" & k==4',
        'label_name in ["N95"] & seq_type=="p" & k==5',
        'label_name in ["N95"] & seq_type=="p" & k==6',
        'label_name in ["N95"] & seq_type=="p" & k==7',
        'label_name in ["N95"] & seq_type=="m1" & k==3',
        'label_name in ["N95"] & seq_type=="m1" & k==4',
        'label_name in ["N95"] & seq_type=="m1" & k==5',
        'label_name in ["N95"] & seq_type=="m1" & k==6',
        'label_name in ["N95"] & seq_type=="m1" & k==7',
        'label_name in ["B05"]',
        'label_name in ["B05"] & seq_type=="p"',
        'label_name in ["B05"] & seq_type=="m1"',
        'label_name in ["B05"] & seq_type=="p" & k==3',
        'label_name in ["B05"] & seq_type=="p" & k==4',
        'label_name in ["B05"] & seq_type=="p" & k==5',
        'label_name in ["B05"] & seq_type=="p" & k==6',
        'label_name in ["B05"] & seq_type=="p" & k==7',
        'label_name in ["B05"] & seq_type=="m1" & k==3',
        'label_name in ["B05"] & seq_type=="m1" & k==4',
        'label_name in ["B05"] & seq_type=="m1" & k==5',
        'label_name in ["B05"] & seq_type=="m1" & k==6',
        'label_name in ["B05"] & seq_type=="m1" & k==7',
    )
    filters_exps_5_bn = (
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95"]',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95"] & seq_type=="p"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95"] & seq_type=="m1"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95"] & seq_type=="p" & k==3',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95"] & seq_type=="p" & k==4',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95"] & seq_type=="p" & k==5',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95"] & seq_type=="p" & k==6',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95"] & seq_type=="p" & k==7',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95"] & seq_type=="m1" & k==3',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95"] & seq_type=="m1" & k==4',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95"] & seq_type=="m1" & k==5',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95"] & seq_type=="m1" & k==6',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["N95"] & seq_type=="m1" & k==7',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["B05"]',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["B05"] & seq_type=="p"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["B05"] & seq_type=="m1"',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["B05"] & seq_type=="p" & k==3',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["B05"] & seq_type=="p" & k==4',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["B05"] & seq_type=="p" & k==5',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["B05"] & seq_type=="p" & k==6',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["B05"] & seq_type=="p" & k==7',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["B05"] & seq_type=="m1" & k==3',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["B05"] & seq_type=="m1" & k==4',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["B05"] & seq_type=="m1" & k==5',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["B05"] & seq_type=="m1" & k==6',
        'algorithm in ["NBK", "BayesNet"] & label_name in ["B05"] & seq_type=="m1" & k==7',
    )
    filters_main_figure = (
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6',
    )
    filters_main_figure_final = (
        'algorithm in ["NBk"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kA" & k==3',
        'algorithm in ["NBk"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kN" & k==5',
        'algorithm in ["BN"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kA" & k==6',
        'algorithm in ["BN"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kN" & k==6',
    )
    filters_main_figure_final_2 = (
        'algorithm in ["NBk"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kA" & k==6',
        'algorithm in ["NBk"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kN" & k==6',
        'algorithm in ["BN"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kA" & k==6',
        'algorithm in ["BN"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kN" & k==6',
    )
    filters_main_figure_final_3 = (
        'algorithm in ["NBk"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kA" & k==3',
        'algorithm in ["NBk"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kN" & k==3',
        'algorithm in ["BN"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kA" & k==3',
        'algorithm in ["BN"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kN" & k==3',
    )
    filters_main_figure_final_4 = (
        'algorithm in ["NBk"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kA"',
        'algorithm in ["NBk"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kN"',
        'algorithm in ["BN"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kA"',
        'algorithm in ["BN"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kN"',
    )
    filters_main_figure_final_s_1 = (
        'algorithm in ["BN"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kA" & k==6',
        'algorithm in ["BN"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70"] & sequence_type=="kN" & k==6',
    )
    filters_final = (
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==3',
        'algorithm in ["NBK"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==5',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="p" & k==6',
        'algorithm in ["BayesNet"] & label_name in ["N95", "N90", "N85", "N80", "N75", "N70"] & seq_type=="m1" & k==6',
    )
    filters_final_2 = (
        'algorithm in ["BN", "DT", "kNN", "SVM"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B05", "B10", "B15", "B20", "B25", "B30"]',
        #'sequence_type == "kN" & k==5 & algorithm in ["BN", "DT", "kNN", "SVM"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B05", "B10", "B15", "B20", "B25", "B30"]',
        #'sequence_type == "kA" & k==3 & algorithm in ["BN", "DT", "kNN", "SVM"] & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B05", "B10", "B15", "B20", "B25", "B30"]',
    )
    filters_final_3 = (
        'model in ["NBk"] & sequence_type=="kA" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B95", "B90", "B85", "B80", "B75", "B70"]',
        'model in ["NBk"] & sequence_type=="kN" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B95", "B90", "B85", "B80", "B75", "B70"]',
        'model in ["tNBk"] & sequence_type=="kA" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B95", "B90", "B85", "B80", "B75", "B70"]',
        'model in ["tNBk"] & sequence_type=="kN" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B95", "B90", "B85", "B80", "B75", "B70"]',
    )
    filters_final_4 = (
      'k==5 & model in ["tNBk"] & sequence_type=="kN" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B95", "B90", "B85", "B80", "B75", "B70"]',
      'k==3 & model in ["tNBk"] & sequence_type=="kA" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B95", "B90", "B85", "B80", "B75", "B70"]',
    )
    filters_final_5 = (
        'model in ["NBk"] & sequence_type=="kN" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B05", "B10", "B15", "B20", "B25", "B30"]',
        'model in ["NBk"] & sequence_type=="kA" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B05", "B10", "B15", "B20", "B25", "B30"]',
        'model in ["tNBk"] & sequence_type=="kN" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B05", "B10", "B15", "B20", "B25", "B30"]',
        'model in ["tNBk"] & sequence_type=="kA" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B05", "B10", "B15", "B20", "B25", "B30"]',
        'k==5 & model in ["tNBk"] & sequence_type=="kN" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B05", "B10", "B15", "B20", "B25", "B30"]',
        'k==3 & model in ["tNBk"] & sequence_type=="kA" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B05", "B10", "B15", "B20", "B25", "B30"]',
        'k==5 & model in ["NBk"] & sequence_type=="kN" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B05", "B10", "B15", "B20", "B25", "B30"]',
        'k==3 & model in ["NBk"] & sequence_type=="kA" & expression_level in ["T95", "T90", "T85", "T80", "T75", "T70", "B05", "B10", "B15", "B20", "B25", "B30"]',
    )
    filters_kdm_1_1 = (
        'all',
        'k==3',
        'k==4',
        'km_set=="10"',
        'km_set=="20"',
        'km_set=="30"',
        'k==3 & km_set=="10"',
        'k==3 & km_set=="20"',
        'k==3 & km_set=="30"',
        'k==4 & km_set=="10"',
        'k==4 & km_set=="20"',
        'k==4 & km_set=="30"',
    )
    filters_kdm_1_2 = (
        'k==4 & expression_type=="RA"',
        'k==4 & expression_type=="PA"',
        'k==3 & expression_type=="RA"',
        'k==3 & expression_type=="PA"',
    )
    filters_kdm_1_3 = (
        'k==3',
        'k==4',
        'k==3 & km_set=="30"',
        'k==3 & expression_type=="PA"',
    )
    filters_ra_1_1 = (
        'all',
        'alphabet=="Bio"',
    )

    classifiers = ('BN_', 'SMO_LG_s', 'DT_s', 'IBk_s', 'BN6_r', 'DT_r', 'IBk_r', 'SMO_LG_r')
    #classifiers = ('BN6', 'SMO_LG', 'DT', 'IBk')
    #classifiers = ('BN6_rmv', 'DT_rmv', 'IBk_rmv', 'SMO_LG_rmv')
    #classifiers = ('SMO_LG_rmv',)
    cur_date = datetime.datetime.now()

    if plot_type == 'n':
        dst_dir = settings.RESULT_DIR + "1_figures/" + cur_date.strftime("%Y-%m-%d")
        #exp_list = exp_list_nbk
        exp_list = exp_list_nbk_1_3
        features = features_1
    elif plot_type == 'n2':
        dst_dir = settings.RESULT_DIR + "1_figures/" + cur_date.strftime("%Y-%m-%d")
        exp_list = exp_list_nbk_final
        features = features_final
    elif plot_type == 'n3':
        dst_dir = settings.RESULT_DIR + "1_figures/" + cur_date.strftime("%Y-%m-%d")
        #exp_list = exp_list_kmd
        #features = features_kmd
        exp_list = exp_list_ra
        features = features_ra
    elif plot_type == 't':
        dst_dir = settings.RESULT_DIR + "1_figures/" + cur_date.strftime("%Y-%m-%d")
        exp_list = exp_list_nbk_1_3
        features = features_1
    elif plot_type == 'r':
        dst_dir = settings.NBK2_RESULT_DIR + "1_figures/" + cur_date.strftime("%Y-%m-%d")
        exp_list = exp_list_rda
        features = features_2
    dir_num = 1

    if os.path.exists(dst_dir):
        dst_dir_new = dst_dir + "-" + str(dir_num)

        while os.path.exists(dst_dir_new):
            dir_num += 1
            dst_dir_new = dst_dir + "-" + str(dir_num)
        dst_dir = dst_dir_new + "/"
    else:
        dst_dir += "/"

    for exp in exp_list:
        if plot_type in ['n', 'n2', 'n3']:
            src_dir = settings.RESULT_DIR + exp + "/1_summary/"
        elif plot_type == 'r':
            src_dir = settings.NBK2_RESULT_DIR + exp + "/1_summary/"
        dst_dir_exp = dst_dir + exp + "/"

        file_list = FileList.ls(path=src_dir, recursion=False, mode=settings.FLS_FILE_ONLY)
        for f_idx, file in enumerate(file_list):
            box_plot = BoxPlot(src=file.full)
            #box_plot.set(showfliers=False)
            #value_list = ['g']
            #box_plot.filter(col_name='label_type', value_list=value_list)
            #box_plot.filter(query='label_type in ["g"]')

            #box_plot.filter(query='seq_type in ["p"]')
            #box_plot.filter(query='seq_type in ["m1"]')
            #box_plot.filter(query='seq_type in ["m1"] & algorithm in ["J48"]')
            #box_plot.filter(query='label_name in ["N95", "N90", "N85", "N80", "N75", "N70"]')

            #for filter_idx, filter_query in enumerate(filters_exps):
            #for filter_idx, filter_query in enumerate(filters_exps_3):
            #for filter_idx, filter_query in enumerate(filters):
            #for filter_idx, filter_query in enumerate(filters_bn):
            #for filter_idx, filter_query in enumerate(filters_exps_4):
            #for filter_idx, filter_query in enumerate(filters_exps_5):
            #for filter_idx, filter_query in enumerate(filters_exps_5_bn):
            #for filter_idx, filter_query in enumerate(filters + filters_exps_5 + filters_bn + filters_exps_5_bn):
            #for filter_idx, filter_query in enumerate(filters_exps_4_1):
            #for filter_idx, filter_query in enumerate(filters_exps_4_2):
            #for filter_idx, filter_query in enumerate(filters + filters_exps_5):
            #for filter_idx, filter_query in enumerate(filters + filters_exps_5 + filters_exps_4):
            #for filter_idx, filter_query in enumerate(filters_nbk + filters_bn):
            #for filter_idx, filter_query in enumerate(filters_main_figure):
            #for filter_idx, filter_query in enumerate(filters_main_figure_final):
            #for filter_idx, filter_query in enumerate(filters_main_figure_final_2):
            #for filter_idx, filter_query in enumerate(filters_main_figure_final_3):
            #for filter_idx, filter_query in enumerate(filters_main_figure_final_4):
            #for filter_idx, filter_query in enumerate(filters_final_2):
            #for filter_idx, filter_query in enumerate(filters_final_3):
            #for filter_idx, filter_query in enumerate(filters_final_4):
            #for filter_idx, filter_query in enumerate(filters_kdm_1_1):
            #for filter_idx, filter_query in enumerate(filters_kdm_1_2):
            #for filter_idx, filter_query in enumerate(filters_kdm_1_3):
            for filter_idx, filter_query in enumerate(filters_ra_1_1):
            #for filter_idx, filter_query in enumerate(filters_final_5):
                print(filter_query)
                if filter_query != 'all':
                    box_plot.filter(query=filter_query)

                for t_idx, feature in enumerate(features):
                    cf_name = 'NBK'
                    for cf in classifiers:
                        if file.name.find(cf) != -1:
                            cf_name = cf
                            break
                    #dst_file = "Figure_{}-{}-{:02d}_{}.png".format(exp, cf_name, t_idx + 1, "_".join(feature))
                    dst_file = "Figure_{}-filter{:02d}-{:02d}_{}.png".format(exp, filter_idx + 1, t_idx + 1, "_".join(feature))
                    #print(dst_file)

                    if not os.path.exists(dst_dir_exp):
                        os.makedirs(dst_dir_exp)
                    #box_plot.show_legend(False)
                    #box_plot.show(features=feature, save_only=True, figure_file=dst_dir_exp + dst_file)
                    box_plot.show(features=feature, save_only=True, figure_file=dst_dir_exp + dst_file, showfliers=False)
                    print(dst_dir_exp + dst_file)


if __name__ == '__main__':
    main(sys.argv)
