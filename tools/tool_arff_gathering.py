"""
    Arff file Gathering Tool
    AUthor: Kyoung Tak Cho
    Created: Sun Apr  5 14:54:56 CDT 2020
    Updated: Tue Apr  8 00:32:48 CDT 2020
"""

import os
import sys
import argparse
import settings
from utils.Common import StrTool, VersionManager
from utils.FileManager import FvManager, FileList, File, FileManager
from shutil import copyfile

parser = argparse.ArgumentParser()
parser.add_argument('-V', '--version', action='version', version=settings.DEV_VERSION,
                    help='Show version.')
parser.add_argument('-s', '--src', help='Source directory path')
parser.add_argument('-d', '--dst', help='Destination directory path. Default is /2_arffs/')
#parser.add_argument('-k', '--ksize', help='K-mer size. If this is set, the original file name will be renamed.')


def main(argv):
    src = None
    dst = None
    k_size = None

    args = parser.parse_args(argv[1:])

    if len(argv) <= 1:
        parser.parse_args(['--help'])
        return

    ## get kmer size
    #if args.ksize:
    #    k_size = args.ksize


    # get source dir
    if args.src:
        src = args.src
        if not os.path.exists(src):
            raise FileNotFoundError('{} does not exist.'.format(src))
        # set target dir
        # get destination dir
        if args.dst:
            dst = args.dst
        else:
            dst = src + '../0_arffs/'      # default

    print('version:', settings.DEV_VERSION)
    print("source directory:", src)
    print("destination directory:", dst)

    #
    # Rename and aggregating to the destination dir
    #
    all_dirs_info = FvManager.get_feature_dirs(src)
    all_dirs_info.sort()
    for dir in all_dirs_info:
        version_info = VersionManager(version_str=dir)
        src_arffs_dir = src + dir + '/0/1/'
        files_list = FileList.ls(path=src_arffs_dir, recursion=False,
                                 mode=settings.FLS_FILE_ONLY, ext=settings.FLS_EXT_ARFF)
        for file in files_list:
            if file.get_ext().lower() == settings.FLS_EXT_ARFF:
                #print(file.full)
                new_file = dst + file.name
                FileManager.dir_create(dst)
                copyfile(src=file.full, dst=new_file)
                print(new_file)


def rename(src, k_size, version, feature_name):
    dst = 'maize_gp_K{KM}_{VER}_{FN}.arff'.format(
        KM=k_size, VER=version, FN=feature_name)

    os.rename(src=src, dst=dst)


if __name__ == "__main__":
    main(sys.argv)
