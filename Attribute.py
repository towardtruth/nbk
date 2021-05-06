'''
    Attribute class
    Author: Kyoung Tak Cho
    Created: Sat Oct 26 07:16:14 CDT 2019
'''


class Attribute(object):
    name = None
    loc = None  # full path
    type = None
    date = None
    owner = None
    description = None

    def set_name(self, name=None):
        if name is None:
            raise ValueError('name is empty.')
        self.name = name

    def set_loc(self, loc=None):
        self.loc = loc

    def set_type(self, type=None):
        self.type = type

    def set_date(self, date=None):
        self.date = date

    def set_owner(self, owner=None):
        self.owner = owner

    def set_description(self, description=None):
        self.description = description


class FileAttribute(Attribute):
    ext = None
    version = None
    contents = None


class FvAttribute(FileAttribute):
    label_type = None
    tissue = None

    def __init__(self, file_name=None, path=None, label_type=None,
                 tissue=None, version=None):
        if file_name is None:
            raise ValueError('file name is empty.')
        self.name = file_name
        self.loc = path
        self.label_type = label_type
        self.tissue = tissue
        self.version = version


class ArffAttribute(FileAttribute):
    #relation = '@relation'
    #attribute = '@attribute'
    #att_class = '@attribute class'
    #data = '@data'
    data = None

    def __init__(self, file_name=None, loc=None, contents=None, data=None):
        if file_name is None:
            raise ValueError('file name is empth.')
        if loc is None:
            raise ValueError('loc is empty.')
        if contents is None:
            raise ValueError('contents is empty.')
        self.name = file_name
        self.loc = loc
        self.contents = contents


