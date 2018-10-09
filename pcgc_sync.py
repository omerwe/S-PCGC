import numpy as np; np.set_printoptions(precision=3, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200)
import os
import time
import sys
import logging
import pcgc_utils

def pcgc_sync(args):

    #read annotations file
    df_annotations = pcgc_utils.load_dfs(args.annot, args.annot_chr, 'annot.gz', 'annot', 'annot', join_axis=0, index_col='SNP')
    df_annotations.drop(columns=['CHR', 'CM', 'BP'], inplace=True)
    
    #read MAFs file
    df_frq = pcgc_utils.load_dfs(args.frqfile, args.frqfile_chr, 'frq', 'frq', 'frqfile', join_axis=0, index_col='SNP', usecols=['SNP', 'MAF'])
    df_frq = df_frq['MAF']
    
    #sync annotations and MAFs file
    if not (df_frq.index.isin(df_annotations.index)).all():
        raise ValueError('Not all SNPs in the annotations file have MAFs')
    is_same = (df_annotations.shape[0] == df_frq.shape[0]) and (df_annotations.index == df_frq.index).all()
    if not is_same:
        df_frq = df_frq.loc[df_annotations.index]
        assert (df_annotations.index == df_frq.index).all()
        
    #create annotations of common SNPs
    is_common = (df_frq>=0.05) & (df_frq <= 0.95)
    df_annotations_common = df_annotations[is_common]
    
    #create overlap matrices
    df_overlap = df_annotations.T.dot(df_annotations)
    df_overlap_common = df_annotations_common.T.dot(df_annotations_common)
    
    #compute min_annot
    min_annot = df_annotations.min(axis=0).values
    min_annot[min_annot>0]=0
        
    #create df_sync
    df_sync = pd.DataFrame(index=df_annotations.columns)
    df_sync['min_annot'] = min_annot
    df_sync['is_continuous'] = [(len(np.unique(df_annotations[c])) > 2) for c in df_annotations.columns]
    df_sync['M'] = df_annotations.sum(axis=0)
    df_sync['M_5_50'] = df_annotations_common.sum(axis=0)
    df_sync['M2'] = np.einsum('ij,ij->j',df_annotations, df_annotations)
    df_sync['M2_5_50'] = np.einsum('ij,ij->j',df_annotations_common, df_annotations_common)
    
    #create fields for non-negative annotations
    df_annotations -= min_annot
    df_annotations_common -= min_annot
    df_sync['M_noneg'] = df_annotations.sum(axis=0)
    df_sync['M_5_50_noneg'] = df_annotations_common.sum(axis=0)
    df_sync['M2_noneg'] = np.einsum('ij,ij->j',df_annotations, df_annotations)
    df_sync['M2_5_50_noneg'] = np.einsum('ij,ij->j',df_annotations_common, df_annotations_common)
    
    df_sync.index.name = 'Category'
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
    
    df_sync, df_overlap, df_overlap_common = pcgc_sync(args)
    df_sync.to_csv(args.out+'.sync', sep='\t', float_format='%0.5e')
    df_overlap.to_csv(args.out+'.overlap', sep='\t', float_format='%0.5e')
    df_overlap_common.to_csv(args.out+'.overlap_5_50', sep='\t', float_format='%0.5e')
    
    
    
        