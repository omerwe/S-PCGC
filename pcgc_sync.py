import numpy as np; np.set_printoptions(precision=3, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200)
import os
import time
import sys
import logging
import pcgc_utils
from tqdm import tqdm




def collect_annotations_chr_info(annot_file, df_frq=None, min_annot=None):
    df_annot = pd.read_table(annot_file, delim_whitespace=True, index_col='SNP', usecols=(lambda c: c not in ['CHR', 'CM', 'BP']))
    
    #compute sync df
    df_sync = pd.DataFrame(index=df_annot.columns)
    df_sync['min_annot'] = np.minimum(df_annot.min(axis=0).values, 0)
    
    #in the first round, we only collect min_annot...
    if min_annot is None:
        return df_sync
        
    #compute overlap matrix and sync fields
    df_overlap = df_annot.T.dot(df_annot)    
    df_sync['is_continuous'] = [(len(np.unique(df_annot[c])) > 2) for c in df_annot.columns]
    df_sync['M'] = df_annot.sum(axis=0)
    df_sync['M2'] = np.einsum('ij,ij->j',df_annot, df_annot)
    df_annot -= min_annot
    df_sync['M_noneg'] = df_annot.sum(axis=0)
    df_sync['M2_noneg'] = np.einsum('ij,ij->j',df_annot, df_annot)
    df_annot += min_annot
    
    #align df_annot and df_frq
    if not (df_annot.index.isin(df_frq.index)).all():
        raise ValueError('Not all SNPs in the annotations file have MAFs')
    is_same = (df_annot.shape[0] == df_frq.shape[0]) and (df_annot.index == df_frq.index).all()
    if not is_same:
        df_frq = df_frq.loc[df_annot.index]
        assert (df_annot.index == df_frq.index).all()
        
    #keep only common SNPs
    is_common = (df_frq>=0.05) & (df_frq <= 0.95)
    df_annot = df_annot[is_common]    
    
    #collect fields for common SNPs
    df_overlap_common = df_annot.T.dot(df_annot)    
    df_sync['M_5_50'] = df_annot.sum(axis=0)
    df_sync['M2_5_50'] = np.einsum('ij,ij->j',df_annot, df_annot)
    df_annot -= min_annot
    df_sync['M_5_50_noneg'] = df_annot.sum(axis=0)
    df_sync['M2_5_50_noneg'] = np.einsum('ij,ij->j',df_annot, df_annot)
    df_annot += min_annot
    
    return df_sync, df_overlap, df_overlap_common


def pcgc_sync2(args):
    
    #define the list of annotation file names
    if args.annot is None and args.annot_chr is None:
        raise ValueError('you must use either --annot or --annot-chr')
    if args.annot is not None and args.annot_chr is not None:
        raise ValueError('you cannot use both --annot and --annot-chr')
    if args.annot is not None:
        annot_fname_list = [args.annot+'annot.gz']
    else:
        annot_fname_list = [args.annot_chr+'%d.annot.gz'%(chr_num) for chr_num in range(1,23)]
    
    #round 1: Collect min_annot fields
    logging.info('Computing annotation minimum values')
    min_annot_list = []
    for annot_fname in tqdm(annot_fname_list, disable=len(annot_fname_list)==1):
        min_annot_list.append(collect_annotations_chr_info(annot_fname).values[:,0])
        
    #compute min_annot across all annotations
    min_annot = np.min(np.array(min_annot_list), axis=0)
    min_annot[min_annot>0] = 0
    
    #collect MAFs
    df_frq = pcgc_utils.load_dfs(args.frqfile, args.frqfile_chr, 'frq', 'frq', 'frqfile', join_axis=0, index_col='SNP', usecols=['SNP', 'MAF'])
    df_frq = df_frq['MAF']
    
    #round 2: collect all fields
    overlap_matrix = 0
    overlap_matrix_common = 0
    df_sync_list = []
    logging.info('Collecting annotation details')
    for annot_fname in tqdm(annot_fname_list, disable=len(annot_fname_list)==1):
        df_sync, df_overlap, df_overlap_common = collect_annotations_chr_info(annot_fname, df_frq, min_annot)
        df_sync_list.append(df_sync)
        overlap_matrix += df_overlap.values
        overlap_matrix_common += df_overlap_common.values
        
    #create df_overlap and df_overlap_common
    df_overlap = pd.DataFrame(overlap_matrix, index=df_overlap.index, columns=df_overlap.columns)
    df_overlap_common = pd.DataFrame(overlap_matrix_common, index=df_overlap_common.index, columns=df_overlap_common.columns)
        
    #Group df_sync results
    func_dict = {}
    func_dict['min_annot'] = np.min
    func_dict['is_continuous'] = np.any
    func_dict['M'] = np.sum
    func_dict['M2'] = np.sum
    func_dict['M_noneg'] = np.sum
    func_dict['M2_noneg'] = np.sum
    func_dict['M_5_50'] = np.sum
    func_dict['M2_5_50'] = np.sum
    func_dict['M_5_50_noneg'] = np.sum
    func_dict['M2_5_50_noneg'] = np.sum    
    assert np.all([c in func_dict for c in df_sync.columns])    
    df_sync_concat = pd.concat(df_sync_list, axis=0)
    df_sync = df_sync_concat.groupby(df_sync_concat.index).agg(func_dict)
    df_sync = df_sync.loc[df_sync_list[0].index]

    #add df index names
    df_sync.index.name = 'Category'
    df_overlap.index.name = 'Category'
    df_overlap_common.index.name = 'Category'
    
    return df_sync, df_overlap, df_overlap_common
    
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--annot', help='prefix of annotations file')
    parser.add_argument('--annot-chr', help='prefix of annotations file in multi-chromosome format')
    parser.add_argument('--frqfile', help='prefix of MAFs file')
    parser.add_argument('--frqfile-chr', help='prefix of MAFs file in multi-chromosome format')
    parser.add_argument('--out', required=True, help='output file prefix')
    args = parser.parse_args()
    
    #check that the output directory exists
    if os.path.isabs(args.out) and not os.path.exists(os.path.dirname(args.out)):    
        raise ValueError('output directory %s doesn''t exist'%(os.path.dirname(args.out)))
        
    pcgc_utils.configure_logger(args.out)
    
    df_sync, df_overlap, df_overlap_common = pcgc_sync2(args)
    df_sync.to_csv(args.out+'.sync', sep='\t')
    df_overlap.to_csv(args.out+'.overlap', sep='\t')
    df_overlap_common.to_csv(args.out+'.overlap_5_50', sep='\t')
    
    
    
        