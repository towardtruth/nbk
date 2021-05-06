"""
    Weka filters for arff
    Author: Kyoung Tak Cho
    Created: Friday, May 15, 2020 9:50:06 PM
    Updated: Friday, May 15, 2020 9:50:06 PM
"""

import os
import sys
import argparse
import settings
from utils.Common import ListTool
from utils.FileManager import FileList, File

parser = argparse.ArgumentParser()
parser.add_argument('-V', '--version', action='version', version=settings.DEV_VERSION,
                    help='Show version.')
parser.add_argument('-s', '--src', help='Source directory')
parser.add_argument('-d', '--dst', help='Destination directory')
parser.add_argument('-f', '--filter', help='Filter type')


def main(argv):
    src = None
    dst = None
    args = parser.parse_args(argv[1:])

    if len(argv) <= 1:
        parser.parse_args(['--help'])
        return

    # get source dir
    if args.src:
        src = args.src
        if not os.path.exists(src):
            raise FileNotFoundError('{} dost not exist.'.format(src))

    # get destination dir
    if args.dst:
        dst = args.dst
    else:
        dst = src + '../0_arffs_rmv/'

    # create dst directory if not exist
    if not os.path.exists(dst):
        os.makedirs(dst)

    #
    # Show information
    #
    print('Version:', settings.DEV_VERSION)
    print('Source Directory:', src)
    print('Destination Directory:', dst)

    #sample: java -classpath "C:\Program Files\Weka-3-8-4\weka.jar" weka.filters.unsupervised.attribute.ReplaceMissingValues -i C:\Users\ktcho\Documents\0_TAKS\research\nbk\5_experiments\0_3_EXP-5\0_arffs\maize_gp_1447_K3_gt01-N95.arff -o C:\Users\ktcho\Documents\0_TAKS\research\nbk\5_experiments\0_3_EXP-5\0_arffs\maize_gp_1447_K3_gt01-N95_rmv.arff

    # get all arffs from source dir
    file_list = FileList.ls(path=src, recursion=False,
                            mode=settings.FLS_FILE_ONLY, ext=settings.FLS_EXT_ARFF)
    cmd_list = list()
    for file in file_list:
        file_src = file.name
        file_src_wo_ext = file.full_wo_ext(file.name)
        file_src_ext = file.get_ext(file.name)
        file_dst = file_src_wo_ext + '_rmv' + file_src_ext

        full_src = src + file_src
        full_dst = dst + file_dst
        #print(full_src)
        #print(full_dst)

        cmd = "java -classpath \"C:\Program Files\Weka-3-8-4\weka.jar\" weka.filters.unsupervised.attribute.ReplaceMissingValues -i \"{}\" -o \"{}\"".format(full_src, full_dst)
        #print(cmd)
        cmd_list.append(cmd)
        #break

    # write batch file
    ListTool.list2file(list_src=cmd_list, file_name=dst + 'ReplaceMissingValues.bat')





if __name__ == "__main__":
    main(sys.argv)


