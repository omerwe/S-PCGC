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
    df_sync = df_annotations.min(axis=0)
    df_sync.loc[df_sync>0]=0
    df_sync.index.name = 'Category'

    return df_sync
    
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--annot', help='prefix of annotations file')
    parser.add_argument('--annot-chr', help='prefix of annotations file in multi-chromosome format')
    parser.add_argument('--out', required=True, help='output file prefix')
    args = parser.parse_args()
    
    #check that the output directory exists
    if os.path.isabs(args.out) and not os.path.exists(os.path.dirname(args.out)):    
        raise ValueError('output directory %s doesn''t exist'%(os.path.dirname(args.out)))
        
    pcgc_utils.configure_logger(args.out)
    
    df_sync = pcgc_sync(args)
    df_sync.to_csv(args.out+'.sync', sep='\t', float_format='%0.5e')
    
    
    
        