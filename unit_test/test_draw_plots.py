"""
    UNIT TEST - Draw Plots
    Author: Kyoung Tak Cho
    Created: Monday, April 13, 2020 9:28:17 AM
    Updated: Monday, April 20, 2020 2:20:10 PM
"""
import unittest
import settings
import os
import datetime
from utils.DataAnalysis import BoxPlot
from utils.FileManager import FileList


class TestDrawPlots(unittest.TestCase):
    def test_draw_box_plot(self):
        src_file = "0_3_EXP-5_BN_10_summary.csv"
        src_path = settings.RESULT_DIR + "0_3_EXP-5/1_summary/" + src_file
        box_plot_exp5_bn = BoxPlot(src=src_path)
        #box_plot_exp5_bn.show(x="K", y="F_measure")
        #box_plot_exp5_bn.show(x="tissue_num", y="F_measure")
        #box_plot_exp5_bn.show(x="label_name", y="F_measure")
        #box_plot_exp5_bn.show(x="label_name", y="F_measure", hue="label_type")
        box_plot_exp5_bn.show(x="tissue_num", y="F_measure", hue="label_type")
        #box_plot_exp5_bn.show(x="tissue_num", hue="label_name", y="F_measure")
        #box_plot_exp5_bn.show(x="tissue_num", hue="label_name", y="F_measure")

    def test_save_box_plot(self):
        '''
        exp_list = ('0_3_EXP-1',
                    '0_3_EXP-2',
                    '0_3_EXP-3',
                    '0_3_EXP-4',
                    '0_3_EXP-5',
                    '0_3_EXP-6',
                    '0_3_EXP-7',
                    '0_3_EXP-8')
        '''
        exp_list = ('0_3_EXP-1-2',)
        features = [('tissue_num', 'f_measure', 'label_type'),
                    ('tissue_num', 'f_measure', 'label_name'),
                    ('tissue_num', 'f_measure', 'k'),
                    ('label_name', 'f_measure', 'k'),
                    ('label_name', 'f_measure', 'label_type')]
        '''
        exp_list = ('0_RA-01',
                    '0_RA-02',)
        #features = [ ('alphabet_size', 'f_measure', 'k')]
        features = [('tissue_num', 'f_measure', 'label_type'),
                    ('alphabet_size', 'f_measure', 'k'),
                    ('alphabet_size', 'f_measure', 'label_type'),
                    ('alphabet_size', 'f_measure', 'label_name'),
                    ('label_name', 'f_measure', 'alphabet_size'),
                    ('alphabet_size', 'f_measure', 'tissue_num'),
                    ('k', 'f_measure', 'alphabet_size'),
                    ('tissue_num', 'f_measure', 'alphabet_size'),
                    ('tissue_num','f_measure', 'label_name'),
                    ('tissue_num', 'f_measure', 'k'),
                    ('label_name', 'f_measure', 'k'),
                    ('label_name', 'f_measure', 'label_type')]
        '''
        classifiers = ('BN_', 'SMO_LG_s', 'DT_s', 'IBk_s', 'BN6_r', 'DT_r', 'IBk_r', 'SMO_LG_r')
        #classifiers = ('BN6', 'SMO_LG', 'DT', 'IBk')
        #classifiers = ('BN6_rmv', 'DT_rmv', 'IBk_rmv', 'SMO_LG_rmv')
        #classifiers = ('SMO_LG_rmv',)
        cur_date = datetime.datetime.now()

        dst_dir = settings.RESULT_DIR + "1_figures/" + cur_date.strftime("%Y-%m-%d")
        #dst_dir = settings.NBK2_RESULT_DIR + "1_figures/" + cur_date.strftime("%Y-%m-%d")
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
            print(exp)
            src_dir = settings.RESULT_DIR + exp + "/1_summary/"
            #src_dir = settings.NBK2_RESULT_DIR + exp + "/1_summary/"
            #dst_dir = src_dir + "../2_figures/"
            dst_dir_exp = dst_dir + exp + "/"

            file_list = FileList.ls(path=src_dir, recursion=False, mode=settings.FLS_FILE_ONLY)
            for f_idx, file in enumerate(file_list):
                box_plot = BoxPlot(src=file.full)
                for t_idx, feature in enumerate(features):
                    cf_name = 'NBK'
                    for cf in classifiers:
                        if file.name.find(cf) != -1:
                            cf_name = cf
                            break
                    dst_file = "Figure_{}-{}-{:02d}_{}.png".format(exp, cf_name, t_idx + 1, "_".join(feature))
                    #print(dst_file)

                    if not os.path.exists(dst_dir_exp):
                        os.makedirs(dst_dir_exp)
                    #box_plot.show(features=feature, save_only=True, figure_file=dst_dir_exp + dst_file, palette='Set3')
                    box_plot.show(features=feature, save_only=True, figure_file=dst_dir_exp + dst_file, palette='husl')
                    print(dst_dir_exp + dst_file)


if __name__ == '__main__':
    unittest.main()
