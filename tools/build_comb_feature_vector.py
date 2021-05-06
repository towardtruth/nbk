'''
    Build Combined Feature Vector
    Author: Kyoung Tak Cho
    Created: Thu Oct 24 01:57:12 CDT 2019
    Updated: Thu Oct 31 11:25:09 CDT 2019
'''

import settings
import sys
import argparse
from utils.Common import StrTool
from utils.FileManager import FvManager

parser = argparse.ArgumentParser()
parser.add_argument('-V', '--version', action='version', version=settings.DEV_VERSION,
                    help='''Show version.''')
parser.add_argument('-c', '--comb', nargs='+',
                    help='''Build Combinations of given features. Example: 1;2,3;4,5;6''')
parser.add_argument('-a', '--all',
                    help='''Build all Combinations of peptide seq and promoter seq''')


def main(argv):
    args = parser.parse_args(argv[1:])

    f_set = list()

    if len(argv) <= 1:
        parser.parse_args(['--help'])
        return

    # get features number to combine
    if args.comb:
        arg_str = args.comb[0]
        f_set = StrTool.str2list_2d(arg_str)
        print('Selected feature set: {}'.format(f_set))

    # generate all combinations of peptide seq and promoter seq
    if args.all:
        pep_builds_not = [392, 393, 394, 398, 399, 400]
        pep_builds = list()
        for i in range(368, 404):
            pep_builds.append(i)
        pep_builds = list(set(pep_builds) - set(pep_builds_not))
        pep_builds.sort()
        #print(pep_builds)

        # build combinations
        for idx, build in enumerate(range(591, 741)):
            f_set.append([pep_builds[idx % len(pep_builds)], build])

        #for target_builds in f_set:
        #    print(target_builds)
        #return

    # build combined feature vector
    loc = settings.RESULT_DIR_COMB
    for target_builds in f_set:
        print('Merging {}'.format(target_builds))
        fv_manager = FvManager(loc=loc, target_builds=target_builds)
        fv_manager.merge()

        # write .arff files
        fv_manager.write_arff()


if __name__ == '__main__':
    main(sys.argv)

