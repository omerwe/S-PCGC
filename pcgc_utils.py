import numpy as np; np.set_printoptions(precision=3, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200)
import os
import time
import sys
import logging
from tqdm import tqdm


def configure_logger(out_prefix):

    logFormatter = logging.Formatter("[%(levelname)s]  %(message)s")
    logger = logging.getLogger()
    logger.setLevel(logging.NOTSET)
    
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)
    
    fileHandler = logging.FileHandler(out_prefix+'.log')
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)


def load_dfs(file_prefix, file_prefix_chr, file_suffix, arg_name, flag_name, join_axis=0, index_col=None, header='infer', use_tqdm=True):

    if file_prefix is not None:
        assert file_prefix_chr is None, 'cannot specify both --%s and --%s-chr'%(flag_name, flag_name)
        logging.info('reading %s file...'%(arg_name))
        fh_list = _splitp(file_prefix)
    elif file_prefix_chr is not None:
        assert file_prefix is None, 'cannot specify both --%s and --%s-chr'%(flag_name, flag_name)
        logging.info('reading %s files...'%(arg_name))
        fh_list = _splitp(file_prefix_chr)
    else:
        raise ValueError('must specify either --%s or --%s-chr'%(flag_name, flag_name))
    if len(fh_list) > 1:
        assert join_axis is not None, 'join_axis cannot be None with multiple file prefixes'
        
    #iterate over file prefixes
    df_outer_list = []
    for fh in fh_list:
        df_fh_list = []
        if file_prefix is not None:
            flist = ['%s%s'%(fh, file_suffix)]
        else:
            flist = ['%s%d.%s'%(fh, chr_num, file_suffix) for chr_num in xrange(1,23)]
            
        #iterate over chromosomes
        df_chr_list = []
        for fname in tqdm(flist, disable=(file_prefix is not None or (not use_tqdm))):
            if not os.path.exists(fname):
                raise IOError('%s not found'%(fname))
            df_chr = pd.read_table(fname, delim_whitespace=True, index_col=index_col, header=header)
            df_chr_list.append(df_chr)
            
        if join_axis is None:
            df = df_chr_list
        else:
            df = pd.concat(df_chr_list, axis=join_axis)
        df_outer_list.append(df)                
            
    #return output
    if len(df_outer_list) == 1: return df_outer_list[0]        
    
    n0 = df_outer_list[0].shape[0]
    for df in df_outer_list[1:]:
        if df.shape[0] != n0:
            raise ValueError('different number of rows found in: %s'%(str(fh_list)))        
    df = pd.concat(df_outer_list, axis=1)
    return df
    
    
def _splitp(fstr):
    flist = fstr.split(',')
    flist = [os.path.expanduser(os.path.expandvars(x)) for x in flist]
    return flist

    
    
    
def find_df_column(df, strings_to_find):
    
    if isinstance(strings_to_find, basestring):
        strings_to_find = [strings_to_find]
        
    is_relevant_col = np.zeros(df.shape[1], dtype=np.bool)
    for str_to_find in strings_to_find:
        is_relevant_col = is_relevant_col | (df.columns.str.upper() == str_to_find.upper())
    if np.sum(is_relevant_col)==0:
        raise ValueError('No matching column found among: %s'%str(strings_to_find))
    elif np.sum(is_relevant_col)>1:
        raise ValueError('Too many matching columns found among: %s'%str(strings_to_find))
    else:
        return df.columns[is_relevant_col][0]
            
        