'''
    Class assignment for building dynamic dataset
    file: build_datasets.py
    Author: Kyoung Tak Cho
    Created: Mon Jul 23 23:25:05 CDT 2018
    Updated: Mon Jul 23 23:25:05 CDT 2018
'''

import settings
import utils.DBManager as dbm
import sqls
from itertools import chain

class ClassAssign(object):
    """ ClassAssign
        Assign class in dataset

        Attributes:
        Methods:
            assign_class(sql)
                return [gene_id, class]
    """

    @staticmethod
    def _assign_class(self, sql):
        """ assign class
        :param sql:
            sql: condition(s)
        :return:
            [gene_id, class]
        """
        # get all gene_ids from walley_full_data|express_abundance(Jesse's dataset)

        # get gene_ids by class(es)

        # assign class(es) for all gene_ids

    @staticmethod
    def _get_all_gene_ids(value_type='p'):
        sql = sqls.get_all_gnid
        par_list = [value_type]
        all_gene_ids = dbm.Pgsql.Common.select_data(sql, par_list)
        print (len(all_gene_ids))

    @staticmethod
    def _low_gene_exp(sql=None):
        if sql is None:
            sql = sqls.sql_low_exp_gene_ids
        low_gene_exp_gene_ids = dbm.Pgsql.Common.select_data(sql, 'g')
        return low_gene_exp_gene_ids

    @staticmethod
    def _class_assign_high_low_exp(gene_prot='g', class_code='l', class_num=2, sql=None, description=None):
        # get gene ids
        gene_ids_raw = dbm.Pgsql.Common.select_data(sql)
        gene_ids_raw.sort()    # sort by gnid
        gene_ids = list(chain.from_iterable(gene_ids_raw))

        class_code = "{GP}{CC}".format(GP=gene_prot, CC=class_code)
        desc = description + sql
        assigned_genes = ",".join(map(str, gene_ids))

        # insert
        dbm.Pgsql.Common.insert_data(sqls.class_assign_new,
                                     (class_code, class_num, desc, assigned_genes))

    @staticmethod
    def gene_low_exp():     # TODO: need to think about value_type option
        # set arguments
        gene_prot = 'g'
        class_code = 'l'
        class_num = 2
        sql = sqls.get_low_exp_gene_ids % gene_prot
        description = "name: low expressed genes,\n\nsql: "
        # new class assignment
        ClassAssign._class_assign_high_low_exp(gene_prot=gene_prot, class_code=class_code, class_num=class_num,
                                               sql=sql, description=description)
    @staticmethod
    def gene_high_exp():
        # set arguments
        gene_prot = 'g'
        class_code = 'h5'
        class_num = 2
        sql = sqls.get_high_exp_gene_ids % (gene_prot, gene_prot)
        description = "name: high expressed genes (top 5%),\n\nsql: "
        # new class assignment
        ClassAssign._class_assign_high_low_exp(gene_prot=gene_prot, class_code=class_code, class_num=class_num,
                                               sql=sql, description=description)
    @staticmethod
    def gene_high_exp_t10():
        # set arguments
        gene_prot = 'g'
        class_code = 'h10'
        class_num = 2
        sql = sqls.get_high_exp_gene_ids_t10 % (gene_prot, gene_prot)
        description = "name: high expressed genes (top 10%),\n\nsql: "
        # new class assignment
        ClassAssign._class_assign_high_low_exp(gene_prot=gene_prot, class_code=class_code, class_num=class_num,
                                               sql=sql, description=description)
    @staticmethod
    def gene_tissues():
        # set arguments
        gene_prot = 'g'
        class_num = 2
        for tissue in range(1,24):
            class_code = 't%s' % tissue
            sql = sqls.get_gene_tissues % (gene_prot, tissue)
            description = "name: genes by tissue# {TS} (top 10%),\n\nsql: ".format(TS=tissue)
            ClassAssign._class_assign_high_low_exp(gene_prot=gene_prot, class_code=class_code, class_num=class_num,
                                                   sql=sql, description=description)



