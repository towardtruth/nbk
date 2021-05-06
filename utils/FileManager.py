'''
    File Manager
    Author: Kyoung Tak Cho
    Created: Fri Jun 15 23:04:43 CDT 2018
    Updated: Sat Nov  2 10:50:24 CDT 2019
'''

from Bio import SeqIO
import csv
import os
import re
import errno
import settings
from utils.Common import VersionManager, DictTool, ListTool
from Attribute import ArffAttribute


class FileManager(object):

    # def __init__(self):
    #    self.conn = None
    #    self.cur = None

    def csv2list(self, filename):
        with open(filename, newline='') as f:
            reader = csv.reader(f)
            rows = []
            for row in reader:
                cols = []
                for col in row:
                    cols.append(col)
                rows.append(cols)
        return rows

    def list2csv(self, listdata, filename):
        f = open(filename, 'w+')
        for row in listdata:
            Line = ",".join(str(value) for value in row)
            f.write("%s\n" % Line)
        f.close()

    def dict2csv(self, dictdata, filename):
        f = open(filename, 'w+')
        for s in dictdata:
            f.write('%s,%d\n' % (s, dictdata[s]))
        f.close()

    @staticmethod
    def file_open(file_name, mode):
        if not os.path.exists(os.path.dirname(file_name)):
            try:
                os.makedirs(os.path.dirname(file_name))
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
        # open file descriptor
        f = open(file_name, mode)
        return f

    @staticmethod
    def dir_create(dir_path):
        if not os.path.exists(dir_path):
            try:
                os.makedirs(dir_path)
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise

    @staticmethod
    def remove_ext(file_name):
        if file_name is None:
            return None
        pos = file_name.rfind('.')
        return file_name[:pos]


#
# Get files/dirs for given path
#
class FileList(object):

    @staticmethod
    def ls(path=None, recursion=False, mode=settings.FLS_ALL, ext=settings.FLS_EXT_ALL):
        if path is None:
            raise ValueError('path is empty.')
        res_list = list()
        for root, dirs, files in os.walk(path):
            if mode == settings.FLS_DIR_ONLY or mode == settings.FLS_ALL:
                for name in dirs:
                    cFile = File(name=name, root=root, f_type=settings.FD_DIR)
                    res_list.append(cFile)
            if mode == settings.FLS_FILE_ONLY or mode == settings.FLS_ALL:
                for name in files:
                    cFile = File(name=name, root=root, f_type=settings.FD_FILE)
                    if ext is not settings.FLS_EXT_ALL:
                        if cFile.is_ext(ext):
                            res_list.append(cFile)
                    else:
                        res_list.append(cFile)
            if not recursion:
                break
        return res_list


#
# File
#
class File(object):
    name = None
    root = None
    full = None
    ext = None
    wo_ext = None
    f_type = None         # either settings.FD_DIR or settings.FD_FILE

    def __init__(self, name=None, root=None, f_type=None):
        self.set_info(name, root, f_type)

    def set_info(self, name=None, root=None, f_type=None):
        if name is not None:
            self.name = name
        if root is not None:
            self.root = root
        if f_type is not None:
            self.type = f_type
        self.full = self.root + self.name

    def print(self):
        for field in self.__dict__:
            print(field)

    def _parse_detail_info(self):
        matched_object1 = re.search(settings.RE_FEATURE_VECTOR_DETAIL_INFO, self.name)
        if matched_object1:
            self.tissue = matched_object1.group()[3:5]
            self.label_type = matched_object1.group()[1:2]
            self.label_name = matched_object1.group()[6:-2]

    def get_ext(self, path=None):
        self._split_ext(path)
        return self.ext.lower()

    def _split_ext(self, path=None):
        if path is None:
            path = self.full
        file_name, file_ext = os.path.splitext(path)
        self.wo_ext = file_name
        self.ext = file_ext

    def full_wo_ext(self, path=None):
        self._split_ext(path)
        return self.wo_ext

    def is_ext(self, ext=None):
        if ext is None:
            return -1

        if self.get_ext() == ext:
            return True
        else:
            return False

    def fs_details(self):
        return File.fs_detail_info(string=self.name)

    def tissue(self):
        return File.tissue_info(string=self.name)

    def version(self):
        return File.version_info(string=self.name)

    @staticmethod
    def tissue_info(string=None):
        return File._search(pattern=settings.RE_SUMMARY_DETAIL_INFO, src_str=string, start=2, cnt=2)

    @staticmethod
    def version_info(string=None):
        return File._search(pattern=settings.RE_FEATURE_VECTOR_VERSION, src_str=string, start=1)

    @staticmethod
    def fs_detail_info(string=None):
        return File._search(pattern=settings.RE_FS_DETAIL_INFO, src_str=string, start=0)

    @staticmethod
    def label_type(string=None):
        return File._search(pattern=settings.RE_SUMMARY_DETAIL_INFO, src_str=string, start=0, cnt=1)

    @staticmethod
    def label_name(string=None):
        return File._search(pattern=settings.RE_SUMMARY_DETAIL_INFO, src_str=string, start=5)

    #@staticmethod
    #def kmer_size(string=None):
    #    return File._search(pattern=settings.RE_K_SIZE, src_str=string, start=5)

    @staticmethod
    def _search(pattern=None, src_str=None, start=None, end=None, cnt=None):
        if pattern is None:
            raise ValueError('pattern is empty.')
        if src_str is None:
            raise ValueError('src_str is empty.')

        if cnt is not None:
            end = start + cnt
        res = None
        matched_object = re.search(pattern=pattern, string=src_str)
        #print(len(matched_object.group()))
        if matched_object:
            res = matched_object.group()[start:end]
        return res


#
# Feature vector file handler
#
class FvFile(object):
    file_name = None
    path = None
    label_type = None
    label_name = None
    tissue = None
    version = None
    kmer_size = None
    data = list()

    def __init__(self, file_name=None, path=None, label_type=None,
                 label_name=None, tissue=None, version=None, kmer_size=None, data=None):
        self.set_info(file_name, path, label_type, label_name, tissue, version, kmer_size, data)
        self._parse_detail_info()

    def set_info(self, file_name=None, path=None, label_type=None,
                 label_name=None, tissue=None, version=None, kmer_size=None, data=None):
        if file_name is not None:
            self.file_name = file_name
        if path is not None:
            self.path = path
        if label_type is not None:
            self.label_type = label_type
        if label_name is not None:
            self.label_name = label_name
        if tissue is not None:
            self.tissue = tissue
        if version is not None:
            self.version = version
        if kmer_size is not None:
            self.kmer_size = kmer_size
        if data is not None:
            self.data = data

    def _parse_detail_info(self):
        #matched_object1 = re.search('_[gpb]t[0-9][0-9][-].*[_][0-9]', self.file_name)
        matched_object1 = re.search(settings.RE_FEATURE_VECTOR_DETAIL_INFO, self.file_name)
        #matched_object2 = re.search('_[0-9].*[.].*[.].*[.].*[0-9]', self.file_name)
        matched_object2 = re.search(settings.RE_FEATURE_VECTOR_VERSION, self.file_name)
        if matched_object1:
            self.tissue = matched_object1.group()[3:5]
            self.label_type = matched_object1.group()[1:2]
            self.label_name = matched_object1.group()[6:-2]
        if matched_object2:
            self.version = matched_object2.group()[1:]

    def print(self):
        print('file name: {},\tpath: {},\tlabel_type: {},\tlabel_name: {},\ttissue: {},\tversion: {}'.format(
            self.file_name, self.path, self.label_type,
            self.label_name, self.tissue, self.version))


class FvFileSet(object):
    def __init__(self):
        self.files = list()  # FvFile * FEATURE_SIZE

        # init feature vector file set
        self._init_set()

    def __iter__(self):
        for file in self.files:
            yield file

    def __len__(self):
        return len(self.files)

    def _init_set(self):
        for i in range(0, settings.FEATURE_SIZE):
            self.files.append(i)

    def add(self, idx=None, fv_file=None):
        if idx is None:
            raise ValueError('idx is empty.')
        if fv_file is None:
            raise ValueError('fv_file is empty.')
        if not isinstance(fv_file, FvFile):
            raise TypeError('fv_file is not an instance of FvFile')
        if idx >= settings.FEATURE_SIZE or idx < 0:
            raise ValueError('idx is out of range.')
        self.files[idx] = fv_file

    def get_file(self, idx=None):
        if idx is None:
            raise ValueError('idx is empty.')
        return self.files[idx]

    def print(self):
        for fv_file in self.files:
            fv_file.print()


class FvManager(object):
    def __init__(self, loc=None, target_builds=None):
        if loc is None:
            raise ValueError('Base location is empty.')
        if target_builds is None:
            raise ValueError('Target builds is empty.')

        self.loc = loc
        self.target_builds = target_builds
        self.fv_dirs = list()
        self.fv = list()  # FvFileSet container : FvFileSet * # of target builds
        self.fv_merged = list()

        # init
        self._get_fv_dirs()
        self._build_file_set()

    def _get_fv_dirs(self):
        if self.target_builds is None:
            raise ValueError('target build is empty.')

        fv_dirs = list()
        feature_dirs = self.find_feature_dirs()
        for build in self.target_builds:
            # check existence of given builds
            for candidate in feature_dirs:
                c_build = self.get_build_num(dir=candidate)
                if c_build == build:
                    fv_dirs.append(self.loc + candidate)
        self.fv_dirs = fv_dirs
        return fv_dirs

    def _build_file_set(self):
        sub_path = '/0/1/'
        #flag_fv = '^feature_vector'
        flag_fv = settings.RE_FEATURE_VECTOR_START

        for fv_dir in self.fv_dirs:
            file_set = FvFileSet()
            final_dir = fv_dir + sub_path
            children = os.walk(final_dir)

            for child in children:
                child[2].sort()
                # ret_files = [file_name for file_name in childs[2] if re.search(flag_fv, file_name)]
                for file_name in child[2]:
                    if re.search(flag_fv, file_name):
                        new_file = FvFile(file_name=file_name, path=final_dir)
                        tissue = new_file.tissue
                        file_set.add(idx=int(tissue)-1, fv_file=new_file)
                break
            self.fv.append(file_set)
            del file_set

    #
    # Merge Feature vectors
    #   return: merged feature vector for 23 tissues as a list() type
    #
    def merge(self):
        fv_tissue = list()

        # build merged fv list by tissue: fv_tissue
        for fv_set in self.fv:
            for idx, file in enumerate(fv_set):
                if len(fv_tissue) <= idx:
                    fv_tissue.append([file])
                else:
                    fv_tissue[idx].append(file)

        # merge feature vectors' contents
        for files in fv_tissue:     # for all 23 tissues
            merged = list()
            tissue = None
            label_name = ''

            # step 1: set fv contents with list() * dict()
            fv_dicts = list()   # list() * dict() * len(files)
            for file in files:
                file_path = file.path + file.file_name
                tmp_dict = CSVLoader.csv2dict(filename=file_path)
                fv_dicts.append(tmp_dict)
                if tissue is None:
                    tissue = int(file.tissue)
                build_num = VersionManager(file.version).build
                label_name += '-' + file.label_type + file.label_name + '_' + str(build_num)

            # step 2: merge
            # 2.1 get common gnids
            common_gnids = list()
            for fv_dict in fv_dicts:
                gnids = list(fv_dict)
                if len(common_gnids) == 0:
                    common_gnids = gnids.copy()
                else:
                    common_gnids = ListTool.common_items(common_gnids, gnids)
            # 2.2 merge
            for gnid in common_gnids:
                fv_m = list()
                for fv_dict in fv_dicts:
                    fv_m += fv_dict[gnid][:-1]
                new_row = [str(gnid)] + fv_m + ['1' if sum([int(fv[gnid][-1]) for fv in fv_dicts]) >= 2 else '0']
                #new_row = [str(gnid)] + fv_m + ['1' if sum([int(fv[gnid][-1]) for fv in fv_dicts]) >= 1 else '0']
                merged.append(new_row)

            # show total number of gnids and possitive class
            print('total gnids: {}\tpositive class#: {}'.format(len(merged), sum([int(x[-1]) for x in merged])))

            # set FvFile with merged FV
            file_name = 'maize_gp_mt%02d%s_%s.arff' % (tissue, label_name, settings.DEV_VERSION)
            path = settings.RESULT_DIR + settings.DEV_VERSION + '/0/1/'
            self.fv_merged.append(FvFile(file_name=file_name, path=path, data=merged))
            print('\t\t===> ' + path + file_name)

    def write_arff(self):

        for fv_file in self.fv_merged:
            cols_size = len(fv_file.data[0]) - 2
            file_name = fv_file.path + fv_file.file_name
            f = FileManager.file_open(file_name, 'w')

            # arff header
            f.write('@relation maize-gp-%s\n' % FileManager.remove_ext(fv_file.file_name))
            for tissue_num in range(0, cols_size):
                f.write('@attribute t%02d numeric\n' % (tissue_num + 1))
            f.write('@attribute class {1,0}\n')
            f.write('@data\n')
            # write data
            for line in fv_file.data:
                f.write('%s\n' % ','.join(line[1:]))    # remove gnid (1st col)
            f.close()

    def find_feature_dirs(self):
        return FvManager.get_feature_dirs(loc=self.loc)

    @staticmethod
    def get_feature_dirs(loc):
        dirs = os.listdir(loc)
        feature_dirs = [dir for dir in dirs if VersionManager.validate_pattern(dir)]
        return feature_dirs

    def get_build_num(self, dir=None, is_ret_int=True):
        if dir is None:
            raise ValueError('dir is empty.')
        version_info = VersionManager(dir)
        if is_ret_int:
            build_num = int(version_info.build)
        else:
            build_num = version_info.build
        return build_num


class CSVLoader(object):
    @staticmethod
    def csv2list(file_name, skip_header=False, numeric=False):
        with open(file_name, newline='') as f:
            reader = csv.reader(f)
            rows = []
            for row in reader:
                cols = []
                for col in row:
                    if numeric:
                        try:
                            cols.append(float(col))
                        except ValueError:
                            cols.append(0.)
                    else:
                        cols.append(col)
                rows.append(cols)
            if skip_header:
                rows.remove(rows[0])
        return rows

    @staticmethod
    def list2csv(list_data, filename):
        #f = open(filename, 'w+')
        f = FileManager.file_open(filename, 'w')
        for row in list_data:
            Line = ",".join(str(value) for value in row)
            f.write("%s\n" % Line)
        f.close()

    @staticmethod
    def dict2csv(dict_data, filename):
        #f = open(filename, 'w+')
        f = FileManager.file_open(filename, 'w')
        for s in dict_data:
            f.write('%s,%d\n' % (s, dict_data[s]))
        f.close()

    @staticmethod
    def csv2dict(filename, skip_header=False, numeric=False):
        with open(filename, newline='') as f:
            reader = csv.reader(f)
            rows = dict()
            for row in reader:
                cols = list()
                rows[int(row[0])] = row[1:]
        return rows


class FileType(object):
    def get_type(self, file_name):
        file_ext_csv = '.CSV'
        file_ext_arff = '.ARFF'
        file_len = len(file_name)

        pos = file_name.upper().find(file_ext_csv)
        if pos == file_len - len(file_ext_csv):
            return file_ext_csv

        pos = file_name.upper().find(file_ext_arff)
        if pos == file_len - len(file_ext_arff):
            return file_ext_arff

    @staticmethod
    def is_csv(file_name):
        ext = file_name[-3:]
        if ext.upper() == 'CSV':
            return True
        return False

    @staticmethod
    def is_arff(file_name):
        ext = file_name[-4:]
        if ext.upper() == 'ARFF':
            return True
        return False


class ArffManager(ArffAttribute):
    @staticmethod
    def write(file_name=None, d_name=None, attributes=None, data=None):
        if file_name is None:
            raise ValueError('file name is empty')
        if d_name is None:
            raise ValueError('dataset name is empty')
        if attributes is None:
            raise ValueError('attributes are empty')
        if data is None:
            raise ValueError('data is empty.')

        f = FileManager.file_open(file_name, 'w')

        # arff header
        f.write('@relation %s\n' % d_name)
        for attr in attributes:
            f.write('@attribute %s numeric\n' % attr)
        f.write('@attribute class {1,0}\n')
        f.write('@data\n')

        # write data
        for line in data:
            f.write('%s\n' % ','.join(map(str, line[1:])))  # remove gnid (1st col)
        f.close()






