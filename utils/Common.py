'''
    Common utils
        list partitioning
    file: utils/Common.py
    Author: Kyoung Tak Cho
    Created: Wed Sep  5 12:14:38 CDT 2018
    Updated: Thu Jul  4 23:22:07 CDT 2019
'''

import os
import re
import errno
import json
import random
import settings


class ListTool(object):

    @staticmethod
    def _partition(src_list=None, capacity=None):
        if src_list is None:
            error_mesg = 'source list is empty.'
            raise ValueError(error_mesg)
        if capacity is None:
            error_mesg = 'partition size is empty.'
            raise ValueError(error_mesg)

        for i in range(0, len(src_list), capacity):
            yield src_list[i:i + capacity]

    @staticmethod
    def partition_by_group_num(src_list=None, size=2, is_sorted=False):
        #capacity = int(len(src_list) / size)
        #return ListTool._partition(src_list, capacity)
        group = list()
        for i in range(0, size):
            group.append([])
        for idx, value in enumerate(src_list):
            group[idx % size].append(value)
        if is_sorted:
            for i in range(0, size):
                group[i] = sorted(group[i])
        return group

    @staticmethod
    def partition_by_group_capacity(src_list=None, capacity=None):
        return ListTool._partition(src_list, capacity)

    @staticmethod
    def list2csv(list_src, file_name, mode='w'):
        #f = open(file_name, 'w')
        f = FileManager.file_open(file_name, mode)
        for row in list_src:
            line = ",".join(str(value) for value in row)
            f.write("%s\n" % line)
        f.close()

    @staticmethod
    def list2file(list_src, file_name, mode='w'):
        f = FileManager.file_open(file_name, mode)
        for row in list_src:
            f.write("%s\n" % str(row))
        f.close()

    @staticmethod
    def list2str(src_list, sep=','):
        #if sep is None:
        #    res_str = " ".join(map(str, src_list))
        #else:
        #    res_str = sep.join(map(str, src_list))
        #return res_str
        return sep.join(map(str, src_list))

    @staticmethod
    def list2range(src_list):
        if len(src_list) != 2:
            raise ValueError('input list should have len() == 2.')
        return StrTool.str2range(','.join(str(value + 1) if idx == 1 else str(value) for idx, value in enumerate(src_list)))

    @staticmethod
    def str2list(src_str, sep=','):
        #if sep is None:
        #    res_list = list(map(int, src_str.split(' ')))
        #else:
        #    res_list = list(map(int, src_str.split(sep)))
        #return res_list
        return list(map(int, src_str.split(sep)))

    @staticmethod
    def twoD2oneD(list_src=None):
        res = list()
        for row in list_src:
            for col in row:
                res.append(col)
        return res

    @staticmethod
    def remove_duplicates(src=None):
        res = None
        if src is not None:
            res = src.copy()
            res = list(dict.fromkeys(res))  # remove duplicates using dict()
        return res

    @staticmethod
    def common_items(list1, list2):
        return [item for item in list1 if item in list2]

    @staticmethod
    def add_rm_dup(list1, list2):
        res = ListTool.remove_duplicates(list1 + list2)
        return res

    @staticmethod
    def sub_f(list1, list2):
        return [item for item in list1 if item not in list2]

    @staticmethod
    def sub(list1, list2):
        return list(set(list1) - set(list2))

    @staticmethod
    def shuffle(list):
        if list is None:
            raise ValueError('list is empty.')
        res = random.shuffle(list)
        return res


class DictTool(object):
    @staticmethod
    def print(dict=None):
        if dict is None:
            raise ValueError('dict is empty.')
        for key, value in dict.items():
            print(key, value)

    @staticmethod
    def dict2str(dict=None):
        if dict is None:
            raise ValueError('dict is empty.')
        return json.dumps(dict)

    @staticmethod
    def str2dict(str=None):
        if str is None:
            raise ValueError('str is empty.')
        return json.loads(str)


class CsvTool(object):

    @staticmethod
    def dict2csv(dict, file_name):
        #f = open(file_name, 'w')
        f = FileManager.file_open(file_name, 'w')
        for row in dict:
            #line = ",".join(str(col) for col in row)
            #f.write("%s\n" % line)
            f.write('%s,%d\n' % (row, dict[row]))
        f.close()


class FileManager(object):

    @staticmethod
    def file_open(file_name, mode):

        # check path (directory) existence and create one(s) if not exists
        if not os.path.exists(os.path.dirname(file_name)):
            try:
                os.makedirs(os.path.dirname(file_name))
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
        # open file descriptor
        f = open(file_name, mode)
        return f


class StrTool(object):

    @staticmethod
    def str2range(range_str):
        range_list = range_str.split(',')
        if len(range_list) < 1 or len(range_list) > 3:
            raise ValueError('range string is not valid.')

        # cast str to int
        range_list = [int(x) for x in range_list]

        if len(range_list) == 1:
            rng = range(range_list[0])
        elif len(range_list) == 2:
            rng = range(range_list[0], range_list[1])
        elif len(range_list) == 3:
            rng = range(range_list[0], range_list[1], range_list[2])

        return rng

    @staticmethod
    def str2list(src_str, sep=','):
        return list(map(int, src_str.split(sep)))

    @staticmethod
    def str2list_2d(src_str, sep1=',', sep2=';'):
        list1 = src_str.split(sep1)
        list_res = [StrTool.str2list(src_str=comb, sep=sep2) for comb in list1]
        return list_res

    @staticmethod
    def ranges2list(src_str):
        list_ranges = list()
        ranges = StrTool.str2list_2d(src_str=src_str, sep1=',', sep2='-')
        for range_list in ranges:
            list_ranges = ListTool.add_rm_dup(list_ranges, list(ListTool.list2range(range_list)))
        return sorted(list_ranges)


class VersionManager(object):
    def __init__(self, version_str):
        self.version = version_str
        self.major = None
        self.minor = None
        self.revision = None
        self.build = None

        # set version information
        self._set_version()

    def _set_version(self):
        if not VersionManager.validate_pattern(self.version):
            raise ValueError('Invalid version!')

        details = VersionManager.detail_info(self.version)
        self.major = int(details[0])
        self.minor = int(details[1])
        self.revision = int(details[2])
        self.build = int(details[3])

    def build_increase(self):
        if self.build is not None:
            self.build += 1
        self.update_version()

    def update_version(self):
        self.version = "".format('{MJ}.{MN}.{RV}.{BD}',
                                 self.major, self.minor, self.revision, self.build)

    def print(self):
        print(VersionManager.detail_info(version=self.version, return_type=2))

    def get_version(self):
        return self.version

    @staticmethod
    def validate_pattern(version=None):
        if version is None:
            raise ValueError('version is empty.')

        # check pattern with RegEx - 1.2.3.4
        if re.search('^' + settings.RE_VERSION, version):
            return True
        else:
            return False

    @staticmethod
    def detail_info(version=None, return_type=None):
        if not VersionManager.validate_pattern(version=version):
            raise ValueError('Invalid Version number format.')

        details = version.split('.')
        if return_type is None or return_type == 1:     # default - list type
            return details
        elif return_type == 2:                          # dictionary type
            return {'major': details[0], 'minor': details[1], 'revision': details[2], 'build': details[3]}
        elif return_type == 3:                          # string type
            return format('')
        else:
            return details                              # default






