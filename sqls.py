'''
    SQLs
    Author: Kyoung Tak Cho
    Created: Sun Aug  5 23:05:35 CDT 2018
'''

###############################################################################
#
# Class assign
#
###############################################################################

#
# build gene_ids table while walley_data and gene_sequence tables are creating.
#
build_gene_ids = """
  INSERT INTO gene_ids
    (gene_id)
  VALUES ('%s')
  ON CONFLICT (gene_id)
  DO NOTHING"""

#
# build walley data
#
build_walley_data = """
  INSERT INTO walley_data
    (scid, gnid, value_type, ea_value)
  VALUES (%s, (SELECT gnid FROM gene_ids WHERE gene_id = %s), %s, %s)"""

#
# build gene_sequence table for peptide sequence
#
build_gene_sequence_pep = """
  INSERT INTO gene_sequence
    (gnid, gene_sub_id, seq, seq_len, seq_type)
  VALUES ((SELECT gnid FROM gene_ids WHERE gene_id = %s),%s,%s,%s,'p')"""

#
# build gene_sequence table for dna sequence
#
build_gene_sequence_dna = """
  INSERT INTO gene_sequence
    (gnid, seq, seq_len, seq_type)
  VALUES ((SELECT gnid FROM gene_ids WHERE gene_id = %s),%s,%s,'d')"""

#
# build gene_sequence table for promoter sequence (DNA)
#
build_gene_sequence_pmt = """
  INSERT INTO gene_sequence
    (gnid, gene_sub_id, seq, seq_len, seq_type)
  VALUES ((SELECT gnid FROM gene_ids WHERE gene_id = %s),%s,%s,%s,%s)"""

#
# buid kmers table
#
build_kmers_get_seq = '''
  SELECT gsid, seq
  FROM gene_sequence
  ORDER BY gsid'''

build_kmers_get_promoter_seq = '''
  SELECT gsid, LEFT(seq, 1000)
  FROM gene_sequence
  WHERE
    seq_type = 'm1'
  ORDER BY gsid'''

build_kmers_get_promoter_seq_5k = '''
  SELECT gsid, LEFT(seq, 1000)
  FROM gene_sequence
  WHERE
    seq_type = 'm5'
  ORDER BY gsid'''

build_kmers_put_kmers = '''
  INSERT INTO kmers
    (k, gsid, kmer, freq)
  VALUES (%s,%s,%s,%s)'''

build_kmers_put_kmers_promoter = '''
  INSERT INTO kmers_promoter
    (k, gsid, kmer, freq)
  VALUES (%s,%s,%s,%s)'''

get_wd_all_gnid = """
  SELECT gnid
  FROM walley_data
  WHERE
    value_type = '%s'
  GROUP BY gnid
  ORDER BY gnid"""

get_wd_all_gnid_both = """
  SELECT gnid
  FROM walley_data
  GROUP BY gnid
  ORDER BY gnid"""

# get all gnids per tissue and label_data_type
get_wd_all_gnid_per_tissue = """
  SELECT gnid
  FROM walley_data
  WHERE
    value_type = '%s'
    AND scid = %s
    AND ea_value %s 0
  ORDER BY gnid"""

# get all gnid per tissue for both (GE & PA combination) data
get_wd_all_gnid_per_tissue_both = """
  SELECT gnid
  FROM walley_data
  WHERE
    scid = %s
    AND ea_value %s 0
  GROUP BY gnid
  ORDER BY gnid"""


# get wd all gnid for debug mode
get_wd_all_gnid_debug = """
  SELECT gnid
  FROM walley_data
  WHERE
    value_type = '%s'
  GROUP BY gnid
  ORDER BY gnid LIMIT %s OFFSET %s"""

# get wd all gnid for debug mode for both (GE and PA combination) data
get_wd_all_gnid_debug_both = """
  SELECT gnid
  FROM walley_data
  GROUP BY gnid
  ORDER BY gnid LIMIT %s OFFSET %s"""

# get missing gnids in promoter sequence data
get_missing_gnids_in_promoter = """
  SELECT DISTINCT wd.gnid
  FROM walley_data AS wd
  WHERE
    wd.gnid NOT IN (
      SELECT DISTINCT gnid FROM gene_sequence WHERE seq_type = 'm1'
    )
  ORDER BY wd.gnid
"""

# get feature group information by fgid
get_feature_group_info_by_fgid = """
  SELECT feature_name, class_size, assigned_genes, corresp_tissue
  FROM feature_group
  WHERE
    fgid = %s"""
# get feature information by fid from features table
get_features_info_by_fsid_corresp_tissue = """
  SELECT feature_name, assigned_genes, corresp_tissue, feature_size
  FROM features
  WHERE
    fsid = %s
    AND corresp_tissue = %s"""

###############################################################################

# get low gene expressed genes
get_low_exp_gene_ids = """
  SELECT gnid
  FROM walley_data
  WHERE
    ea_value = 0
    AND value_type = '%s'
  GROUP BY gnid
  HAVING count(*) = 23
  ORDER BY gnid"""

# get high expressed genes (top 5%)
get_high_exp_gene_ids = """
  SELECT gnid
  FROM
    (SELECT gnid, wd.scid, (ea_value / mx.max_value) AS vp
    FROM walley_data wd,
      (SELECT scid, max(ea_value) AS max_value
      FROM walley_data
      WHERE
        value_type = '%s'
      GROUP BY scid) AS mx
    WHERE
      wd.scid = mx.scid
      AND value_type = '%s'
    ORDER BY scid, vp DESC) AS nm
  GROUP BY gnid
  ORDER BY min(nm.vp) DESC
  limit 2189"""

# get high expressed genes (top 10%)
get_high_exp_gene_ids_t10 = """
  SELECT gnid
  FROM
    (SELECT gnid, wd.scid, (ea_value / mx.max_value) AS vp
    FROM walley_data wd,
      (SELECT scid, max(ea_value) AS max_value
      FROM walley_data
      WHERE
        value_type = '%s'
      GROUP BY scid) AS mx
    WHERE
      wd.scid = mx.scid
      AND value_type = '%s'
    ORDER BY scid, vp DESC) AS nm
  GROUP BY gnid
  ORDER BY min(nm.vp) DESC
  limit 4378"""

# get genes by each tissues, top 10%
get_gene_tissues = """
  SELECT gnid
  FROM walley_data
  WHERE
    value_type = '%s'
    AND scid = %s
    AND ea_value > %f
  ORDER BY gnid"""

# get genes by each tissues at random - FOR TEST
get_gene_tissues_random = """
  SELECT gnid
  FROM (
    SELECT gnid
    FROM walley_data
    WHERE
      value_type = '%s'
      AND scid = %s
      AND gnid NOT IN ( %s )
      AND ea_value > 0
    ORDER BY random() limit %s
  ) as genes_random
  ORDER BY gnid"""

# get genes for building dataset of High gene expression Low protein abundance data
get_gene_hglp = """
  SELECT wd.gnid
  FROM walley_data wd,
    (
      SELECT gnid
      FROM walley_data
      WHERE
        value_type = 'g' 
        AND scid = {SCID} 
        AND ea_value > {GCUT}
      ORDER BY gnid
    ) AS nwd
  WHERE
    wd.scid = {SCID}
    AND wd.gnid = nwd.gnid
    AND wd.value_type = 'p'
    AND wd.ea_value < {PCUT}"""

# get genes for building dataset of gene expression and protein abundance combo - ignore zero value
get_gene_gp_comb = """
  SELECT wd.gnid
  FROM walley_data wd,
    (
      SELECT gnid
      FROM walley_data
      WHERE
        value_type = '{GP1}' 
        AND scid = {SCID} 
        AND ea_value {INQ1} {CUT1}
      ORDER BY gnid
    ) AS nwd
  WHERE
    wd.scid = {SCID}
    AND wd.gnid = nwd.gnid
    AND wd.value_type = '{GP2}'
    AND wd.ea_value {INQ2} {CUT2}
    AND wd.ea_value <> 0"""

# get genes from GE dataset with new cutoffs
get_gene_by_cutoffs = """
  SELECT wd.gnid
  FROM walley_data AS wd,
    (
      SELECT *
      FROM cutoffs
      WHERE
        percentile IN ('{PCNT}')
        AND gp_type = '{GP}'
        AND tsid = {TSID}
    ) AS co
  WHERE
    wd.value_type = co.gp_type
    AND wd.scid = co.tsid
    AND wd.ea_value {INEQ} co.cutoff
    AND wd.ea_value <> 0
  ORDER BY wd.gnid"""

###############################################################################
# Feature set & Features
###############################################################################

#
# add new features - NEW (2019-05-20)
#
new_features = """
  INSERT INTO features
    (fsid, feature_name, description, assigned_genes, corresp_tissue, feature_size)
  VALUES (%s, %s, %s, %s, %s, %s)"""


#
# add new feature set with returning added feature set id (fsid)
#
new_feature_set = """
  INSERT INTO feature_set
    (fs_name, label_data_type, description, class_size)
  VALUES (%s, %s, %s, %s)
  RETURNING fsid"""


#
# add new feature group - old (not used, 2019-05-20)
#
new_feature_group = """
  INSERT INTO feature_group
    (feature_name, class_size, description, assigned_genes, corresp_tissue)
  VALUES (%s, %s, %s, %s, %s)"""
  #ON CONFLICT (feature_name)
  #DO NOTHING"""


#
# Get feature set by fsid
#
get_feature_set = """
  SELECT fs_name, label_data_type, class_size
  FROM feature_set
  WHERE
    fsid = %s"""


#
# Get features by fsid
#
get_features = """
  SELECT fid, feature_name, assigned_genes, corresp_tissue, feature_size
  FROM features
  WHERE
    fsid = %s
  ORDER BY corresp_tissue"""

#
# Get Tissue IDs by fsid
#
get_crsp_tissues_fsid = """
  SELECT corresp_tissue
  FROM features
  WHERE
    fsid = %s
  ORDER BY corresp_tissue"""


###############################################################################

# get low protein abundance genes
sql_low_prot_abd_gene_ids = """
  SELECT gnid
  FROM walley_data
  WHERE
    ea_value = 0
    AND value_type = 'p'
  GROUP BY gnid
  HAVING count(*) = 23"""



###############################################################################
# Sequence
###############################################################################

get_gnid_by_gsid = """
  SELECT gnid
  FROM gene_sequence
  WHERE
    gsid = %s"""


get_seq_info_by_gnid_seq_len = """
  SELECT gsid, gene_sub_id, seq
  FROM gene_sequence
  WHERE
    gnid = %s
    AND seq_type = '%s'
  ORDER BY seq_len {MM}
  LIMIT 1"""


###############################################################################
# KmerFreq class
###############################################################################

# get k-mer frequency
get_kmer_freq_by_gnid = """
  SELECT kmer, freq, gsid
  FROM kmers
  WHERE
    gsid = (
      SELECT gsid
      FROM gene_sequence
      WHERE
        gnid = %s
        AND seq_type = '%s'
      ORDER BY seq_len DESC
      LIMIT 1
    )
    AND k = %s"""

get_kmer_freq_by_gsid = """
  SELECT kmer, freq
  FROM %s
  WHERE
    k = %s
    AND gsid = %s
  ORDER BY kmer"""

get_kmer_freq_by_gsid_seq_type = """
  SELECT kmer, freq
  FROM %s
  WHERE
    k = %s
    AND gsid = %s
    AND seq_type = %s
  ORDER BY kmer"""

#
# Get sum of kmer freq for given set of gsids
#
get_kmer_freq_sum = """
SELECT kmer, sum(freq)
FROM kmers
WHERE
    k = %s
    AND gsid in (%s)
GROUP BY kmer
ORDER BY sum(freq) DESC"""

get_kmer_freq_promoter_sum = """
SELECT kmer, sum(freq)
FROM kmers_promoter
WHERE
    k = %s
    AND gsid in (%s)
GROUP BY kmer
ORDER BY sum(freq) DESC"""


###############################################################################
# Cutoffs class
###############################################################################

# get cutoffs for certain percentiles
get_cutoffs_by_percentiles = """
  SELECT gp_type, tsid, percentile, cutoff
  FROM cutoffs
  WHERE
    percentile IN ({PCTS}, 0)
  GROUP BY gp_type, tsid, percentile, cutoff
  ORDER BY tsid, gp_type, percentile"""



###############################################################################
# Dataset table
###############################################################################

build_dataset = """
  INSERT INTO dataset
    (tissue, expression_type, expression_level, name, 
     positive_class, negative_class)
  VALUES (%s, %s, %s, %s, %s, %s)"""

get_dataset = """
  SELECT dsid, tissue, expression_type, expression_level, name, positive_class, negative_class
  FROM dataset
  WHERE
    expression_type in (%s)
    AND expression_level in (%s)
    AND tissue in (%s)
  ORDER BY expression_type, expression_level DESC, tissue
"""

get_seq_by_gsids = """
  SELECT gsid, gnid, seq
  FROM gene_sequence
  WHERE
    gsid in (%s)"""


###############################################################################
# kmer_freq table
###############################################################################

build_kmer_freq = """
  INSERT INTO kmer_freq
    (dsid, class, percentile, freq, kmer, k)
  VALUES (%s, %s, %s, %s, %s, %s)"""

get_kmer_by_top_rank = """
  SELECT kmer
  FROM kmer_freq
  WHERE
    dsid = %s
    AND k = %s
    AND class = '%s'
  ORDER BY freq DESC
  LIMIT %s"""



