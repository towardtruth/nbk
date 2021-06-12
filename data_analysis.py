"""
    Data Analysis

    Author: Kyoung Tak Cho
    Created: Thu Jun 27 12:14:21 CDT 2019
    Updated: Sun Dec  1 22:11:35 CST 2019
"""

import sys
from utils.DataAnalysis import PrecisionRecallCurve, DistributionOfGenes
import settings


def main(argv):

    # draw PR curve
    pr_curve()

    # draw genes distribution graph
    #gene_dist()


def gene_dist():
    plot_genes_dist = DistributionOfGenes()
    plot_genes_dist.load('../data_analysis/distributions/distribution_test_GE_tissue_1.csv')
    plot_genes_dist.plot_gene_dist()


def pr_curve():
    #gp_type = 'g'
    gp_type = 'p'
    tissue_sub_name = '-N70'
    for tissue in range(1, 24):
        pr_curve = PrecisionRecallCurve(tissue)
        #pr_curve.load('../results/1.3.9.167/0/1/prediction_details_pt{}_1.3.9.167.csv'.format(tissue))
        #pr_curve.load('../results/1.3.10.168/0/1/prediction_details_gt{}_1.3.10.168.csv'.format(tissue))
        pr_curve.load('../results/{VER}/0/1/prediction_details_{GP}t{TS:02d}{TSN}_{VER}.csv'.format(
            VER=settings.DEV_VERSION, GP=gp_type, TS=tissue, TSN=tissue_sub_name))
        pr_curve.pr_curve()
        del pr_curve


#
# Run main program
#
if __name__ == '__main__':
    main(sys.argv)
