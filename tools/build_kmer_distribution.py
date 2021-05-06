"""
    Build K-mer Distribution
    Author: Kyoung Tak Cho
    Created: Mon Jan 27 21:09:08 CST 2020
    Updated: Fri Jan 31 09:59:52 CST 2020
"""
import sys
from GeneGroups import Features, FeatureSet
import sqls
from utils.DBManager import Pgsql
from utils.Common import ListTool


def main(argv):
    build_kmer_dist()


def get_missing_gnids_in_promoter():
    sql = sqls.get_missing_gnids_in_promoter
    res = Pgsql.Common.select_data(sql=sql)
    missing_gnids_in_promoter = Pgsql.Common.to_list(res)
    return missing_gnids_in_promoter


def remove_missing_gnids(all_gnids, tissue_id):
    # remove missing gnids for promoter
    missing_gnids = get_missing_gnids_in_promoter()
    # new_sub_dataset = list(set(all_gnids) - set(missing_gnids))
    new_sub_dataset = ListTool.sub(all_gnids, missing_gnids)
    print('Tissue#: {}, all gnids: {}, common items: {}, removed: {}'.format(
        tissue_id,
        len(all_gnids),
        len(ListTool.common_items(all_gnids, missing_gnids)),
        len(new_sub_dataset)))
    return sorted(new_sub_dataset)


def build_kmer_dist():
    use_max_seq_len = True
    seq_tyeps = ['m1', 'p']
    file_name_seq_type = ""
    #kmer_size = 3
    kmer_sizes = [2,3,4,5,6,7]
    fsids = [44,45,46,47,48,49, 50,51,52,53,54,55, 68,69,70,71,72,73, 74,75,76,77,78,79, 83,84,85,86,87,88]

    for seq_type in seq_tyeps:
        for kmer_size in kmer_sizes:
            for fsid in fsids:
                fs_info = FeatureSet.get(fsid)
                fs_name = fs_info.name.strip()
                #seq_type = fs_info.label_data_type.strip()

                kmer_dist = list()
                # add hader
                kmer_dist.append(['FSID', 'Dataset Name', 'Tissue#', 'Kmer', 'Frequency'])

                tids = FeatureSet.tissues(fsid)
                conn = Pgsql.Common.connect()
                for tissue_id in tids:
                    assigned_gnids = Features.get_gnids_pos(fsid, tissue_id, cur=conn.cursor())
                    if seq_type == 'm1':
                        assigned_gnids = remove_missing_gnids(assigned_gnids, tissue_id)

                    # get gsids
                    assigned_gsids = Features.get_gsids_pos(assigned_gnids, seq_type, use_max_seq_len, cur=conn.cursor())
                    assigned_gsids = sorted(assigned_gsids)
                    str_assigned_gsids = ','.join(str(gsid) for gsid in assigned_gsids)

                    # get kmer distribution by seq_type
                    if seq_type == 'p':
                        res = Pgsql.Common.select_data(sql=sqls.get_kmer_freq_sum,
                                                       pars=(kmer_size, str_assigned_gsids),
                                                       cur=conn.cursor())
                        # set file name
                        file_name_seq_type = "protein_seq"
                    elif seq_type == 'm1':
                        res = Pgsql.Common.select_data(sql=sqls.get_kmer_freq_promoter_sum,
                                                       pars=(kmer_size, str_assigned_gsids),
                                                       cur=conn.cursor())
                        # set file name
                        file_name_seq_type = "promoter_seq"

                    # make list for saving result file
                    for kmer_freq in res:
                        kmer_dist.append([fsid, fs_name, tissue_id] + list(kmer_freq))

                # write kmer_dist file (.csv)
                file_name = 'data_anl/{}_K{}_{}_{}.csv'.format(file_name_seq_type, kmer_size, fsid, fs_name)
                ListTool.list2csv(kmer_dist, file_name)
                print('{} saved.'.format(file_name))

                conn.close()


#
# Run main program
#
if __name__ == '__main__':
    main(sys.argv)


