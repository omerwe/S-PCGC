import numpy as np; np.set_printoptions(precision=3, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200)
import os
import time
import sys
import ldsc_parse as ps
import logging
import ldscore_r2
import pcgc_utils


def compute_r2_prod(bfile, annot_prefix, annot_chr_prefix, sync_prefix, ld_wind_cm=1, chunk_size=50):

    # read bim/snp
    array_snps = ps.PlinkBIMFile(bfile+'.bim')
    snpnames = array_snps.df['SNP']

    #read annotations file
    if annot_prefix is None and annot_chr_prefix is None:
        df_annotations = pd.DataFrame(np.ones(snpnames.shape[0]), index=snpnames, columns=['base'])
    else:
        df_annotations = pcgc_utils.load_dfs(annot_prefix, annot_chr_prefix, 'annot.gz', 'annot', 'annot', join_axis=0, index_col='SNP')
        df_annotations.drop(columns=['CHR', 'CM', 'BP'], inplace=True)    
        
        #apply min_annot correction to ensure no negative values
        category_names = df_annotations.columns
        if sync_prefix is None:
            raise ValueError('--annot must be used together with --sync')        
        df_sync = pd.read_table(sync_prefix+'sync', header=None, delim_whitespace=True, index_col=0, squeeze=True)        
        if df_sync.shape[0] != len(category_names) or not np.all(df_sync.index == category_names):
            raise ValueError('Annotations in sync file do not match those in annotations/prodr2 files')
        min_annot = df_sync.values
        df_annotations -= min_annot
    
    #keep only annotations of SNPs in plink file
    if not snpnames.isin(df_annotations.index).all():
        raise ValueError('not all SNPs have annotations')
    df_annotations = df_annotations.loc[snpnames]
    
    #find #individuals in bfile
    df_fam = pd.read_table(bfile+'.fam', header=None)
    n = df_fam.shape[0]
    
    #read plink file    
    keep_indivs = None
    mafMin = None
    reload(ldscore_r2)
    logging.info('Loading SNP file...')
    geno_array = ldscore_r2.PlinkBEDFile(bfile+'.bed', n, array_snps, keep_snps=None,
        keep_indivs=keep_indivs, mafMin=mafMin)
        
    #compute r2_prod_table
    logging.info('Computing r2 prod...')
    coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]
    block_left = ldscore_r2.getBlockLefts(coords, ld_wind_cm)
    if block_left[len(block_left)-1] == 0:
        error_msg = 'Only a single block selected - this is probably a mistake'
        raise ValueError(error_msg)
    t0 = time.time()
    geno_array._currentSNP = 0
    r2prod_table = geno_array.ldScoreVarBlocks(block_left, chunk_size, annot=df_annotations.values)
    
    df_r2prod_table = pd.DataFrame(r2prod_table, index=df_annotations.columns, columns=df_annotations.columns)
    return df_r2prod_table
        


    
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--bfile', required=True, help='plink file')
    parser.add_argument('--annot', default=None, help='prefix of annotations file')
    parser.add_argument('--annot-chr', default=None, help='Multi-chromosome prefix of annotations file')
    parser.add_argument('--sync', default=None, help='Prefix of pcgc sync file, created by pcgc_sync.py')
    parser.add_argument('--out', required=True, help='output file prefix')
    parser.add_argument('--ld-wind-cm', type=float, default=1.0, help='window size to be used for estimating r2 products in units of centiMorgans (cM).')
    parser.add_argument('--chunk-size',  type=int, default=50, help='chunk size for r2 product calculation')
    args = parser.parse_args()
    
    #check that the output directory exists
    if os.path.isabs(args.out) and not os.path.exists(os.path.dirname(args.out)):    
        raise ValueError('output directory %s doesn''t exist'%(os.path.dirname(args.out)))
        
    pcgc_utils.configure_logger(args.out)        
    
    df_prod_r2 = compute_r2_prod(args.bfile, args.annot, args.annot_chr, args.sync, ld_wind_cm=args.ld_wind_cm, chunk_size=args.chunk_size)
    df_prod_r2.to_csv(args.out+'.prodr2', sep='\t', float_format='%0.5e')
    
    
    
        