'''
    Build aggregated result summary
    Author: Kyoung Tak Cho
    Created: Thu Oct 31 11:25:09 CDT 2019
    Updated: Sunday, June 21, 2020 3:33:33 AM
'''

import settings
import os
import sys
import argparse
from utils.Common import StrTool, VersionManager
from utils.FileManager import FvManager, FileList, CSVLoader, File

parser = argparse.ArgumentParser()
parser.add_argument('-V', '--version', action='version', version=settings.DEV_VERSION,
                    help='''Show version.''')
parser.add_argument('-b', '--builds', nargs='+',
                    help='''Build Combinations of given features. Example: 1;2,3;4,5;6''')
#parser.add_argument('-a', '--all',
#                    help='''Build all Combinations of peptide seq and promoter seq''')
parser.add_argument('-s', '--src', help='''source directory path''')
parser.add_argument('-q', '--seq_type', help='''seq type''')


def main(argv):
    builds_list = None
    source_dir_path = None
    main_exp_num = None
    seq_type = None

    args = parser.parse_args(argv[1:])

    if len(argv) <= 1:
        parser.parse_args(['--help'])
        return

    # get seq type
    if args.seq_type:
        seq_type_str = args.seq_type
        if seq_type_str.lower() in ['p', 'protein', 'peptide']:
            seq_type = 'Protein'
        elif seq_type_str.lower() in ['d', 'dna']:
            seq_type = 'DNA'
        elif seq_type_str.lower() in ['m', 'promoter', 'pmt', 'pm']:
            seq_type = 'Promoter'
        elif seq_type_str.lower() in ['r', 'rda', 'ra']:
            seq_type = 'Reduced Alphabet'
        else:
            seq_type = 'Unknown - ' + seq_type_str
        print('seq type: {}'.format(seq_type))

    # get features number to combine
    if args.builds:
        arg_str = args.builds[0]
        print(arg_str)
        if arg_str.lower() == 'all':
            builds_list = None
        else:
            builds_list = StrTool.ranges2list(arg_str)
    else:
        builds_list = None      # default: all
        print('Builds: {}'.format(builds_list))

    # get source directory
    if args.src:
        source_dir_path = args.src
        main_exp_num = source_dir_path[-20:-11]
        if not os.path.exists(source_dir_path):
            raise FileNotFoundError('{} does not exist.'.format(source_dir_path))
    print('source directory:', source_dir_path)

    #
    # build aggregated results summary
    #
    aggregated_summary = None
    exp_num = None
    all_dirs_info = FvManager.get_feature_dirs(source_dir_path)
    all_dirs_info.sort()
    for dir in all_dirs_info:
        version_info = VersionManager(version_str=dir)
        if builds_list is None or version_info.build in builds_list:
            prediction_summary_dir = source_dir_path + dir + '/0/'
            # get prediction summary file
            files_list = FileList.ls(path=prediction_summary_dir, recursion=False, mode=settings.FLS_FILE_ONLY)
            for file in files_list:
                print(file.full)
                #print(file.version())
                # read summary file and get detail information
                #   tissue#, label_type, label_name
                summary_list = CSVLoader.csv2list(file_name=file.full)

                # add overall rows
                for row in summary_list[1:]:
                    if row[3] == 'overall':
                        tissue = File.tissue_info(row[2])
                        label_type = File.label_type(row[2])
                        label_name = File.label_name(row[2])
                        exp_num = version_info.build

                        # set header with additional columns (seq_type, label_type, label_name, tissue#)
                        if aggregated_summary is None:
                            header = summary_list[0]
                            if seq_type == "Reduced Alphabet":
                                header[0:1] = ['Exp_Num', 'Mapping ID']
                                header[3:1] = ['seq_type', 'label_type', 'label_name', 'tissue_num']
                            else:
                                header[0] = 'Exp_Num'
                                header[2:1] = ['seq_type', 'label_type', 'label_name', 'tissue_num']
                            aggregated_summary = list()
                            aggregated_summary.append(header)
                        data = row
                        if seq_type == "Reduced Alphabet":
                            data[0:1] = [exp_num, data[0]]
                            data[3:1] = [seq_type, label_type, label_name, tissue]
                        else:
                            data[0] = exp_num
                            data[2:1] = [seq_type, label_type, label_name, tissue]
                        aggregated_summary.append(data)

    # write aggregated summary
    if seq_type == "Reduced Alphabet":
        file_name = os.path.join(source_dir_path + '../1_summary/',
                             'prediction_summary_aggregated_RA_' + settings.DEV_VERSION + '.csv')
    else:
        file_name = os.path.join(source_dir_path + '../1_summary/',
                             'prediction_summary_aggregated_' + main_exp_num + '_' + settings.DEV_VERSION + '.csv')
    CSVLoader.list2csv(list_data=aggregated_summary, filename=file_name)
    print('{} has been created.'.format(file_name))


if __name__ == "__main__":
    main(sys.argv)


