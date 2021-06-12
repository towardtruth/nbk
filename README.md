# NBK (Naive Bayes k-mer model)
Predicting tissue specific mRNA and protein abundance in maize: A machine learning approach.

## Files & Direcoties structures

```
nbk/
├── data
│   ├── gene.out.csv
│   └── prot.out.csv
├── unit_test
│   └── kmer_test.py
├── utils
│   ├── Common.py
│   ├── DataAnalysis.py
│   ├── DBManager.py
│   └── FileManager.py
├── build_datasets.py
├── build_data_walley_full.py                       NOT USED
├── ClassAssign.py                                  NOT USED
├── data_analysis.py
├── Data.py
├── GeneGroups.py
├── gene_prediction.py
├── KmerFreq.py                                     NOT USED (-> models.py)
├── Kmer.py                                         NOT USED (-> models.py)
├── models.py
├── settings.py
├── sqls.py
└── Validation.py

```

### gene_prediction.py
```
Main program for predicting mRNA and protein abundance (gene expression)
```


### build_datasets.py

```
    Dataset builder (Application)
    Descriptions:
        Build dataset and insert into the database
    
    Possible options:
        --dataset
        --db_type (NOT USED)
        --get_sample_code (NOT COMPLTED YET)
    Datasets:
        Peptide sequence data
        DNA sequence data
        Promotor sequence data
        Walley data (Labeling data)
        Kmers
    
    Usage:
      1. build promoter data
        python build_datasets.py --dataset=pm
      2. build kmers data
        python build_datasets.py --dataset=km
    
    Functions:
        main(argv)
        build_kmers()
        build_pep_seq_data()
        build_dna_seq_data()
        build_pmt_seq_data()
        build_walley_data()
        get_sample_code()
```

### data_anlysis.py

```
    Data Analysis (Application)
    
    Descriptions:
      Tools for Data analysis such as drawing distributions of data and drawing precision recall curve.
    
    Daya anlysis:
      - Distribution of genes 
      - Precision recall curve
    
    Actual utils:
      utils.DataAnalysis
        PrecisionRecallCurve
        DistributionOfGenes
    
    Usages:
      python data_analysis.py
    
    TODO:
      - adding parse to support options in prompt command
    
```


### Data.py
```
    Data builder of genes and feature vector
    
    class GetData:
      Descriptions:
        Build genes dataset from walley data
      Static Methods:
        wd_all_gnid(gene_prot=None, debug_mode=False)
        wd_all_gnid_per_tissue(exp_setting, scid)
        to_list(src)
      
      TODO:
        Refactorization:
          Change class name: GetData -> GetGenes
        Adding:
          GetGenes
          build negative class set
    
    class FeatureVector:
    
    Methods:
      __init__(self, exp_setting)
      __del__(self)
      _set_genes_info(self)
      _set_wd_all_gnids(self)
      _set_features(self)
      _print_feature_dataset
      pg_db_connect(self)
      get_features(self)
      init_feature_dataset(self, fold_size=None)
      validation_small_genes(self, class_size=None, gene_size=None, kmer-size=None)
      cross_validation(self)
      set_confusion_matrix(self, validation=None, fold_size=None)
      write_genes_info(self)
      write_prediction_results(self)
      write_prediction_summary(self)
      write_feature_vector(self)
      create_arff(self)
      build_feature_vector(self)
      build_test(self)
```

### GeneGroups.py
```
    A gene group can be defined by tissues and abundance cutoffs.
    
    FeatureSet class
      methods:
        __init__(self, name=None, fsid=None, label_data_type=None, description=None,
                 class_size=None, features=None, tissues=None)
        _set(self, name=None, fsid=None, label_data_type=None, description=None,
             class_size=None, features=None, tissues=None)
        get(fsid=None)
        tissues(fsid=None)
        exp_level(self)
        exp_type(self)
        
    Features class
      methods:
        __init__(self, target_features=None, gnid_min_max=None, test_mode=False, 
                 exp_setting=None)
        __iter__(self)
        __next__(self)
        __len__(self)
        get_features(self)
        get_target_features(self)
        get_feature(self, idx)
        _init_features(self)
        get_gnids_pos(fsid=None, tid=None, cur=None)
        get_gsids_pos(gnids=None, seq_type=None, get_max_len=True, cur=None)
        get_gnids_neg(feature_set_id=None, seq_type=None, label_type=None, tissue_id=None, 
                      exclusive_gnids=None, size=None)
        remove_missing_gnids(self, gnids)
        get_exclusive_gnids(self, gnids, missing_gnids)
        get_missing_gnids_in_promoter()
        _set_features(self, test_mode=False)
      staticmethods:
        _low_gene_exp(sql=None)
        _new_feature_group(feature_code=None, class_size=2, sql=None, description=None, 
                           corresp_tissue=None, gene_prot='g')
        _new_features(fsid=None, feature_code=None, sql=None, description=None, 
                      corresp_tissue=None, gene_prot='g')
        gene_low_exp()
        gene_high_exp()
        gene_high_exp_t10()  - for TEST only
        gene_tissues()
        prot_tissues(gene_prot=None, class_size=2, fsid=None)
        class_assign_tissues(gene_prot=None, fsid=None)
        class_assign_hglp(gp_type=None, fsid=None)
        class_assign_by_percentile(fsid=None, gp_type=None, percentile=None, is_top=True)
        class_assign_gp_comb(comb_conf=None, exp_setting=None)
        class_assign_GE_N90(gp_type=None, fsid=None)  - for TEST only
        class_assign_PA_N90(gp_type=None, fsid=None)  - for TEST only
        add_new_feature_set(fs_name, gp_type, description, class_size)
        get_cutoff(gene_prot)
        
    Genes class
        __init__(self, gene_ids=None, seq_type='p', k=3)
        __iter__(self)
        __len__(self)
        _set_genes(self, gene_ids)
        get_gene(self, gnid=None)
    
    FreqDictSet class
        __init__(self)
        _init_fd_fold(self)
        _assert_k(self, k=None)
        append_kmer_dict(self, k=None, freq_dict=None)
        set_fd_fold(self, fd_fold=None, k=None)
        get_fd_fold(self, k=None, fold_idx=None)
        get_fd_class(self, k=None)
        get_fd_complement(self, k=None, fold_idx=None)
    
    Partition class - NOT USED
    
```


### models.py
```
models
  - defined all models

  GeneSequence class
  KmerFreq class
    methods:
      __init__(self, seq=None, k=None, conn=None):
        self._set_kmer_freq()
      get_kmer_freq(self)
      _set_kmer_freq(self)
      
  FreqDict class
    methods:
      __init__(self, k=None)
      __add__(self, other)
      __sub__(self, other)
      set_kmer_freq(self, kmer_freq)
      add_kmer_freq(self, kmer, freq)
      sub_kmer_freq(self, kmer, freq)
      get_kmer_freq_value(self, kmer=None)
      total_freq(self)
      total_freq_float(self)
      
  ValidationDataset class
    Description:
      ValidationDataset handles genes per class and all genes divided into 
      the folds by random partitioning.
    fields:
      gnids                 # all gene ids
      fold_size
      group                 # list() type and each element conresponds to per fold
      confusion_matrix
    
    methods:
      __init__(self, gnids=None, fold_size=None):
        self._random_partition()
      __len__(self):
      __iter__(self):
      _random_partition(self):
      get_fold_size(self):
      set_confusion_matrix(self, cm=None, fold_idx=None):
      get_confusion_matrix(self, fold_idx=None):


  Gene class
  GenePrediction class
  ConfusionMatrix class
  Summary class
  Feature class
  GeneGroup class
  SequenceGroup class
  FeatureInfo class
  Configuration class
  Cutoff class
  Cutoffs class
  CombConf class
```

### settings.py
```
    Manage all environmental variables in whole program wide
```

### sqls.py
```
    Manage all SQL queries used in the program
```


### Validation.py
```
Validation package
  ValidationFreqDict class
    methods:
      __init__(self, datasets=None, k=None, genes=None)
        self._init_freq_dict_set()
      _init_freq_dict_set(self)
      get_fd_class(self, class_num=None, k=None)
      get_fd_fold(self, class_num=None, k=None, fold_idx=None)
      get_train_set(self, class_num=None, k=None, fold_idx=None)
      test_freq_dict_set(self)
      
  CrossValidation class
    methods:
      __init__(self, genes=None, all_gnids=None,
               class_size=None, fold_size=None, kmer_size=None,
               exp_setting=None)
      build_datasets(self, assigned_genes=None, neg_class_mode=1, corresp_tissue=None)
      test_datasets(self)
      validation(self)
      get_fold_size(self)
      
  SmallValidation class - NOT USED
```

## Database structure

### walley_data


### cutoffs


### gene_ids

### gene_sequence


### kmers



### feature_set


### features



### sample_code or tissue_code






# Modules by Functionality

## Gene assignment

We have binary classification and it has two classes which are positive and negative class.
When we have assigned genes means that these genes associate to positive class. 
On the other hand, 


## compute/load k-mer frequency
```
Get k-mer frequency data from DB
  class: KmerFreq (models.py)
    _set_kmer_freq()
      1. get kmer & freq from DB
        sqls.get_kmer_freq_by_gsid
        pars: k, seq.get_gsid()
      2. set frequency data into the self.kmer_freq (FreqDict)
        self.kmer_freq.add_kmer_freq(kmer, freq)
Set k-mer frequency data into the kmer_freq using FreqDict
  class: FreqDict
    self.kmer_freq : dict()
    add_kmer_freq(kmer, freq)
      self.kmer_freq[kmer] = self.kmer_freq.get(kmer, 0) + freq
      # update(add) frequency of kmer with freq
```

## Genes instance

```
Genes
  self.genes = dict()
  self.seq_type
  
  _set_genes(gene_ids)
    self.genes[gnid] = gene     # gene : Gene() : GeneSequence() * 3
    
```



# TEST modules

## Unit TEST


## system TEST


### system_test.py


```
Description:


current process: gene_prediction.py 
  - set exp_settings
    set_target_features()
    set_class_size
    set_kmer_size
  - FeatureVector()
    TEST MODE
      test_features_dataset()
    NORMAL MODE
      cross_validation()
      write_prediction_results()
      build_feature_vector()
      write_feature_vector()
      create_arff()
      write_prediction_summary()

Design:
  - select gnids
  - 

```















