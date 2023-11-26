import numpy as np; np.set_printoptions(precision=3, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200)
import scipy.stats as stats
import os
import time
import sys
import logging
from scipy.stats import t as tdist
import string
import random
import pcgc_utils
import scipy.optimize as optimize
from functools import reduce

MIN_NUM_SNPS = 100
    
def splash_screen():
    print('*********************************************************************')
    print('* S-PCGC for heritability and genetic correlation estimates')
    print('* Version 2.0.0')
    print('* (C) 2018 Omer Weissbrod')
    print('*********************************************************************')
    print()
    

class SPCGC_Data:
    def __init__(self, args, sumstats_prefix, sumstats_prefix_chr, category_names, sum_l2):
        
        #load summary statistics
        df_sumstats, var_t, pve, mean_Q, N = self.read_sumstats(args, sumstats_prefix, sumstats_prefix_chr)        

        #read Gty files and compute trace/deflation ratios
        if args.he:
            df_Gty = self.create_synthetic_Gty(N, 1.0, category_names)
            trace_ratios = np.ones(1)
            deflation_ratio = 1.0
        else:            
            df_Gty, M_base = self.read_Gty_files(args, sumstats_prefix, sumstats_prefix_chr, category_names, N, mean_Q)
            mean_l2 = sum_l2 / df_sumstats.shape[0]
            trace_ratios, deflation_ratio = \
                self.compute_deflation_ratios(sumstats_prefix, sumstats_prefix_chr, pve, N, mean_l2, M_base)
        
        #save class members
        self.df_Gty = df_Gty
        self.df_sumstats = df_sumstats
        self.trace_ratios = trace_ratios
        self.deflation_ratio = deflation_ratio
        self.var_t = var_t
        self.mean_Q = mean_Q
        self.N = N
        
        
    def read_Gty_files(self, args, sumstats_prefix, sumstats_prefix_chr, category_names, N, mean_Q):
    
        #read otherstats files
        df_otherstats_list = pcgc_utils.load_dfs(sumstats_prefix, sumstats_prefix_chr, 'otherstats', 'otherstats', 'otherstats', join_axis=None, use_tqdm=False, allow_duplicates=True)

        #sum M2 (sum of squares of annotations) of each annotation
        M_annot_sumstats2 = np.zeros(len(category_names))
        for df_otherstats_chr in df_otherstats_list:            
            for c_i, c in enumerate(category_names):
                if len(category_names)==1:
                    df_M_annot_sumstats2 = df_otherstats_chr.loc[df_otherstats_chr['Property'].str.startswith('M2_'), 'Value']
                else:
                    df_M_annot_sumstats2 = df_otherstats_chr.query('Property == "M2_%s"'%(c))['Value']
                if df_M_annot_sumstats2.shape[0] == 0:
                    raise ValueError('M2_%s not found in otherstats file'%(c))
                if df_M_annot_sumstats2.shape[0] > 1:
                    raise ValueError('Multiple M2_%s values found in otherstats file'%(c))
                M_annot_sumstats2[c_i] += df_M_annot_sumstats2.values[0]

                
        #load PCGC Gty files and aggregate them (across chromosomes)
        if args.no_Gty:                    
            df_Gty = self.create_synthetic_Gty(N, mean_Q, category_names)
        else:
            df_Gty = pcgc_utils.load_dfs(sumstats_prefix, sumstats_prefix_chr, 'Gty.gz', 'Gty', 'Gty', join_axis=0, index_col=['fid', 'iid'], use_tqdm=False, allow_duplicates=True)
            df_Gty = np.sqrt((df_Gty**2).groupby(['fid', 'iid']).sum())
            
            #synchronize df_Gty columns to match category_names
            if len(category_names) > 1:
                columns_intersect = df_Gty.columns.intersection(category_names)
                if len(columns_intersect) < len(category_names):
                    raise ValueError('Gty files and prodr2 files must have the same annotations')
                if len(columns_intersect) < df_Gty.shape[1]:
                    logging.warning('Gty file has ununsed annotations')
                    df_Gty.loc[:, category_names]
                if not np.all(df_Gty.columns == category_names):
                    df_Gty = df_Gty.loc[:, category_names]
                
            #normalize Gty by the total number of SNPs in the genome
            for anno_i in range(df_Gty.shape[1]):
                df_Gty.iloc[:,anno_i] /= np.sqrt(M_annot_sumstats2[anno_i])                
              
        M_base = M_annot_sumstats2[0]
        return df_Gty, M_base
            
                    
    def create_synthetic_Gty(self, N, mean_Q, category_names, str_len=10):
        random_str = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(str_len))
        iid = [random_str+str(i) for i in range(1,int(N+1))]
        df_Gty = pd.DataFrame(iid, columns=['fid'])
        df_Gty['iid'] = iid            
        df_Gty.set_index(['fid','iid'], inplace=True, drop=True)
        for c in category_names:
            df_Gty[c] = np.sqrt(mean_Q)
        return df_Gty
    
    
    def read_sumstats(self, args, sumstats_prefix, sumstats_prefix_chr):

        #load summary statistics
        if args.he:        
            try:
                df_sumstats = pcgc_utils.load_dfs(sumstats_prefix, sumstats_prefix_chr, '', 'sumstats', 'sumstats', index_col='SNP', allow_duplicates=True)
            except IOError:            
                df_sumstats = pcgc_utils.load_dfs(sumstats_prefix, sumstats_prefix_chr, 'sumstats.gz', 'sumstats', 'sumstats', index_col='SNP', allow_duplicates=True)
        else:
            df_sumstats = pcgc_utils.load_dfs(sumstats_prefix, sumstats_prefix_chr, 'sumstats.gz', 'sumstats', 'sumstats', index_col='SNP', allow_duplicates=True)
        
        #transform Z column if it wasn't created especially for PCGC
        if args.he:
            if 'pcgc_sumstat' not in df_sumstats.columns:
                if 'Z' not in df_sumstats.columns:
                    raise ValueError('cannot find a Z column in summary statistics file')
                else:
                    df_sumstats['pcgc_sumstat'] = df_sumstats['Z'] * np.sqrt(df_sumstats['N'])
        
        #if HE regrssion is used, create default naive values for PCGC-relevant fields
        if args.he:
            var_t = 0
            mean_Q = 1
            pve = np.array([])
            N = df_sumstats['N'].mean()
        #load PCGC other stats files
        else:            
            df_otherstats_list = pcgc_utils.load_dfs(sumstats_prefix, sumstats_prefix_chr, 'otherstats', 'otherstats', 'otherstats', join_axis=None, use_tqdm=False, allow_duplicates=True)
            df_otherstats = df_otherstats_list[0]
            var_t = df_otherstats.query('Property == "var_t"')['Value'].values[0]
            pve = df_otherstats.loc[df_otherstats['Property'].str.startswith('pve'), 'Value'].values
            mean_Q = df_otherstats.query('Property == "mean_Q"')['Value'].values[0]
            N = df_otherstats.query('Property == "N"')['Value'].values[0]
        
        #filter out SNPs with very large summary statistics
        if args.chisq_max is None: chisq_max = max(0.001*df_sumstats.N.max(), 80)            
        else: chisq_max = args.chisq_max
        df_sumstats['Z2'] = df_sumstats['pcgc_sumstat']**2 / df_sumstats['N'] / mean_Q
        is_large_z = (df_sumstats['Z2'] > chisq_max)
        if is_large_z.any():
            logging.warning('Removing %d summary statistics with Z^2 > %0.1f'%(is_large_z.sum(), chisq_max))
            df_sumstats = df_sumstats.loc[~is_large_z]
            
        return df_sumstats, var_t, pve, mean_Q, N
    
        
        
    def compute_deflation_ratios(self, sumstats_prefix, sumstats_prefix_chr, pve, N, mean_l2, M_base):
                
        #read diagGRM files
        df_diagGRM = pcgc_utils.load_dfs(sumstats_prefix, sumstats_prefix_chr, 'diagGRM.gz', 'diagGRM', 'diagGRM', join_axis=0, index_col=['fid', 'iid'], use_tqdm=False)
        df_diagGRM = df_diagGRM.groupby(['fid', 'iid']).sum()                
                
        #compute trace ratio
        deflate_columns   = df_diagGRM.columns.str.startswith('diag_G_deflate')
        nodeflate_columns = df_diagGRM.columns.str.startswith('diag_G_nodeflate')
        sum_diag_deflate = df_diagGRM.loc[:,deflate_columns].sum(axis=0).values
        sum_diag_nodeflate = df_diagGRM.loc[:,nodeflate_columns].sum(axis=0).values
        trace_ratios = sum_diag_nodeflate / sum_diag_deflate
        
        #compute deflation ratios        
        if (len(pve) == 0):
            if np.isclose(trace_ratios[0], 1):
                deflation_ratio = 1.0
            else:
                raise ValueError('trace deflation reported, but no pve values found!')
        else:
            var_diag_deflate = df_diagGRM.loc[:,deflate_columns].iloc[:,0].var(ddof=0) / M_base**2
            var_diag_nodeflate = df_diagGRM.loc[:,nodeflate_columns].iloc[:,0].var(ddof=0) / M_base**2
            trace_G = sum_diag_nodeflate[0] / M_base
            deflation_ratio = self.compute_PCS_deflation_ratio(pve, N, M_base, mean_l2,
                                                           trace_G,
                                                           var_diag_nodeflate,
                                                           var_diag_deflate)
        return trace_ratios, deflation_ratio
        
        
    def compute_PCS_deflation_ratio(self, pve, N, M_base, mean_l2,
                               trace_G,                               
                               var_diag_deflate,
                               var_diag_nodeflate):

        #transform pve to eigenvalues
        s = pve * trace_G
        
        #approximate sum of squares of off-diagonals elements of G, using the derivation of Bulik-Sullivan
        sum_offdiag_G2 = N**2 / float(M_base) * mean_l2
        
        #approximate sum of squares of diagonal elements of G, np.sum(np.diag(G))
        sum_diag_G2 = N*var_diag_nodeflate + trace_G**2 / N
        
        #approximate sum of squares of eigenvalues (which is also the sum of squares of G, np.sum(G**2))
        sum_s2 = sum_offdiag_G2 + sum_diag_G2
        
        #approximate sum of squares of Gr (i.e., np.sum(Gr**2))
        sum_Gr2 = sum_s2 - s.dot(s)
        
        #approximate sum of diagonal elements of Gr (i.e., np.sum(np.diag(Gr**2)))
        sum_diag_Gr2 = (trace_G - s.sum())**2 / N   +   N * var_diag_deflate
        
        #approximate sum of squares of off-diagonals elements of Gr
        sum_offdiag_Gr2 = sum_Gr2 - sum_diag_Gr2
        
        #compute deflation ratio
        deflation_ratio = sum_offdiag_Gr2 / sum_offdiag_G2
        
        return deflation_ratio
        
        

class SPCGC_RG:
    def __init__(self, obj_h2_1, obj_h2_2, obj_cov, M_annot, category_names):

        #handle negative h2 entries
        assert len(category_names) == obj_h2_1.cat.shape[0]
        assert len(category_names) == obj_h2_2.cat.shape[0]
        for cat_name, obj_h2_1_c, obj_h2_2_c in zip(category_names, obj_h2_1.cat, obj_h2_2.cat):
            if obj_h2_1_c<0 or obj_h2_2_c<0:
                logging.warning('Cannot compute rg for %s because of negative h2 estimates'%(cat_name))

        #compute overall genetic correlation and its stderr
        if obj_h2_1.tot>=0 and obj_h2_2.tot>=0:
            rg = obj_cov.tot / np.sqrt(obj_h2_1.tot * obj_h2_2.tot)
        else:
            rg = np.nan
        tot_jk_1 = obj_h2_1.delete_values.dot(M_annot)
        tot_jk_2 = obj_h2_2.delete_values.dot(M_annot)
        tot_jk_cov = obj_cov.delete_values.dot(M_annot)
        i = (tot_jk_1>=0) & (tot_jk_2>=0)
        rg_jk = np.zeros(tot_jk_1.shape) + np.nan
        rg_jk[i] = tot_jk_cov[i] / np.sqrt(tot_jk_1[i]*tot_jk_2[i])
        rg_var = np.var(rg_jk, ddof=0) * len(obj_cov.delete_values-1)
        
        #compute per-annotation genetic correlation
        i = (obj_h2_1.cat>=0) & (obj_h2_2.cat>=0)
        rg_annot = np.zeros(obj_h2_1.cat.shape) + np.nan
        rg_annot[i] = obj_cov.cat[i] / np.sqrt(obj_h2_1.cat[i] * obj_h2_2.cat[i])
        cat_annot_jk_1 = obj_h2_1.delete_values * M_annot
        cat_annot_jk_2 = obj_h2_2.delete_values * M_annot
        cat_annot_jk_cov = obj_cov.delete_values * M_annot
        i = (cat_annot_jk_1>=0) & (cat_annot_jk_2>=0)
        rg_annot_jk = np.zeros(cat_annot_jk_1.shape) + np.nan
        rg_annot_jk[i] = cat_annot_jk_cov[i] / np.sqrt(cat_annot_jk_1[i]*cat_annot_jk_2[i])
        rg_annot_cov = np.cov(rg_annot_jk.T, ddof=0) * len(obj_cov.delete_values-1)
        if rg_annot_jk.shape[1] == 1: rg_annot_cov = np.array([[rg_annot_cov]])
        rg_annot_se = np.sqrt(np.diag(rg_annot_cov))
        
        #create an rg df
        self.df_annot_rg = pd.DataFrame(index=category_names)
        self.df_annot_rg['RG'] = rg_annot
        self.df_annot_rg['SE'] = rg_annot_se
        
        #save class members
        self.rg = rg
        self.rg_var = rg_var
        self.rg_var = rg_var
        self.rg_se = np.sqrt(rg_var)
        self.rg_annot = rg_annot
        self.rg_annot_cov = rg_annot_cov
        self.rg_annot_se = rg_annot_se
                
        
class SPCGC_Cov:
    def __init__(self, coef, delete_values, intercept, delete_intercepts, M_annot, category_names, overlap_matrix, M_tot, var_t1, var_t2, is_continuous):
        self.coef = coef
        self.delete_values = delete_values
        self.intercept = intercept
        self.delete_intercepts = delete_intercepts
        self.compute_enrichment_suffstats(M_annot)
        var_t = np.sqrt(var_t1 * var_t2)
        
        #create output df
        n_blocks = self.delete_values.shape[0]
        self.df_enrichment = self._overlap_output(category_names, n_blocks, overlap_matrix, np.row_stack(M_annot).T, M_tot, print_coefficients=True, is_continuous=is_continuous)
        
        #compute total h2 and save it to a df
        h2_marginal, h2_conditional, h2_marginal_se, h2_conditional_se = \
            self.compute_h2(M_annot, var_t)
        df_dicts = []
        df_dicts.append({'Quantity':'h2_conditional', 'Value':h2_conditional, 'SE':h2_conditional_se})
        df_dicts.append({'Quantity':'h2_marginal', 'Value':h2_marginal, 'SE':h2_marginal_se})
        #df_dicts.append({'Quantity':'Intercept', 'Value':self.intercept, 'SE':self.intercept_se})        
        self.df_h2 = pd.DataFrame(df_dicts, columns=['Quantity', 'Value', 'SE'])
            
        
    def compute_enrichment_suffstats(self, M_annot):
        n_blocks = self.delete_values.shape[0]
        self.coef_cov = np.cov(self.delete_values.T, ddof=0) * (n_blocks-1)
        if len(M_annot)==1:
            self.coef_cov = np.array([[self.coef_cov]])
        self.coef_se = np.sqrt(np.diag(self.coef_cov))
        self.cat = self.coef * M_annot
        self.cat_cov = self.coef_cov * np.outer(M_annot,M_annot)
        self.tot = self.cat.sum()
        self.tot_cov = self.cat_cov.sum()
        self.prop = self.cat / self.tot        
        self.prop_cov = np.cov((self.delete_values*M_annot / np.sum(self.delete_values*M_annot, axis=1)[:,np.newaxis]).T, ddof=0) * (n_blocks-1)
        
        #intercept se
        self.intercept_cov = np.var(self.delete_intercepts, ddof=0) * (n_blocks-1)
        self.intercept_se = np.sqrt(self.intercept_cov)
            
        
    def compute_h2(self, M_annot, var_t):
        h2_cond = self.coef.dot(M_annot)
        h2_marg = h2_cond / (1+var_t)
        h2_cond_se = np.sqrt(np.var(self.delete_values.dot(M_annot), ddof=0) * (len(self.delete_values)-1))
        h2_marg_se = h2_cond_se / (1+var_t)                
        return h2_marg, h2_cond, h2_marg_se, h2_cond_se

        
    def _overlap_output(self, category_names, n_blocks, overlap_matrix, M_annot, M_tot, print_coefficients, is_continuous):
        '''LD Score regression based code for summarizing enrichment of overlapping categories.'''
        n_annot = len(category_names)
        overlap_matrix_prop = np.zeros([n_annot,n_annot])        
        M_annot[np.isclose(M_annot, 0)] = 1e-12 #prevent numerical errors...
        for i in range(n_annot):            
            overlap_matrix_prop[i, :] = overlap_matrix[i, :] / M_annot
            
        prop_hsq_overlap = np.dot(
            overlap_matrix_prop, self.prop.T).reshape((1, n_annot))
        prop_hsq_overlap_var = np.diag(
            np.dot(np.dot(overlap_matrix_prop, self.prop_cov), overlap_matrix_prop.T))
        prop_hsq_overlap_se = np.sqrt(
            np.maximum(0, prop_hsq_overlap_var)).reshape((1, n_annot))
        one_d_convert = lambda x: np.array(x).reshape(np.prod(x.shape))
        prop_M_overlap = M_annot / M_tot
        enrichment = prop_hsq_overlap / prop_M_overlap
        enrichment_se = prop_hsq_overlap_se / prop_M_overlap
        overlap_matrix_diff = np.zeros([n_annot,n_annot])
        for i in range(n_annot):
            if not M_tot == M_annot[0,i]:
                overlap_matrix_diff[i, :] = overlap_matrix[i,:]/M_annot[0,i] - \
                    (M_annot - overlap_matrix[i,:]) / (M_tot-M_annot[0,i])

        diff_est = np.dot(overlap_matrix_diff, self.coef)
        diff_cov = np.dot(np.dot(overlap_matrix_diff, self.coef_cov),overlap_matrix_diff.T)
        diff_se = np.sqrt(np.diag(diff_cov))
        diff_p = np.array([np.nan if diff_se[i]==0 else 2*tdist.sf(abs(diff_est[i]/diff_se[i]), n_blocks) \
            for i in range(n_annot)])
        
        zscore = np.zeros(n_annot)
        is_zero_both =  ((self.coef_se==0) & (self.coef==0))
        is_zero_denom = ((self.coef_se==0) & (~is_zero_both))
        is_valid_zscore = (self.coef_se!=0)
        zscore[is_valid_zscore] = one_d_convert(self.coef[is_valid_zscore]) / one_d_convert(self.coef_se[is_valid_zscore])
        zscore[is_zero_both] = 0.0
        zscore[is_zero_denom] = np.nan
        
        #mark out enrichment-related output for continuous annotations
        prop_M_overlap[0,is_continuous] = np.nan
        prop_hsq_overlap[0,is_continuous] = np.nan
        prop_hsq_overlap_se[0,is_continuous] = np.nan
        enrichment[0,is_continuous] = np.nan
        enrichment_se[0,is_continuous] = np.nan
        diff_p[is_continuous] = np.nan
        
        df = pd.DataFrame({
            'Category': category_names,
            'Prop._SNPs': one_d_convert(prop_M_overlap),
            'Prop._h2': one_d_convert(prop_hsq_overlap),
            'Prop._h2_std_error': one_d_convert(prop_hsq_overlap_se),
            'Enrichment': one_d_convert(enrichment),
            'Enrichment_std_error': one_d_convert(enrichment_se),
            'Enrichment_p':diff_p,
            'Coefficient': one_d_convert(self.coef),
            'Coefficient_std_error': self.coef_se,
            'Coefficient_z-score': zscore
            #'Coefficient_z-score': one_d_convert(self.coef) / one_d_convert(self.coef_se)
        })
        if print_coefficients:
            df = df[['Category', 'Prop._SNPs', 'Prop._h2', 'Prop._h2_std_error',
                    'Enrichment','Enrichment_std_error', 'Enrichment_p',
                     'Coefficient', 'Coefficient_std_error','Coefficient_z-score']]
        else:
            df = df[['Category', 'Prop._SNPs', 'Prop._h2', 'Prop._h2_std_error',
                    'Enrichment','Enrichment_std_error', 'Enrichment_p']]
        return df
        
            

        
            

class SPCGC:
    def __init__(self, args):

        #create prodr2 table
        df_prodr2 = self.load_prodr2(args)
        
        #create the initial index of all SNPs
        index_intersect = self.load_all_snp_indices(args)                
                
        #create and sync objects for all the sumstats files
        sum_l2 = df_prodr2.iloc[0,0]
        pcgc_data_list = self.create_study_objects(args, sum_l2, category_names=df_prodr2.columns)
        index_intersect = self.sync_data_files(args, pcgc_data_list, index_intersect)
        
        #load all the SNP-related data files
        M_annot, df_prodr2, df_annotations_sumstats_noneg, df_sync, df_overlap, df_l2, df_w_ld = \
            self.load_annotations_data(args, df_prodr2, index_intersect)
                                
        #compute h2 and gencov of each pair of studies
        gencov_arr = np.empty((len(pcgc_data_list), len(pcgc_data_list)), dtype=object)
        for i in range(len(pcgc_data_list)):
            oi = pcgc_data_list[i]
            for j in range(i+1):
                oj = pcgc_data_list[j]                
                cov_ij = self.create_cov_obj(args, oi, oj,
                                             df_annotations_sumstats_noneg, df_prodr2, df_sync, df_overlap, M_annot,
                                             df_l2=df_l2, df_w_ld=df_w_ld)
                gencov_arr[i,j] = cov_ij
                gencov_arr[j,i] = cov_ij
                               
                               
        #compute rg
        rg_arr = np.empty((len(pcgc_data_list), len(pcgc_data_list)), dtype=object)
        for i in range(len(pcgc_data_list)):
            for j in range(i+1):
                rg_arr[i,j] = SPCGC_RG(gencov_arr[i,i], gencov_arr[j,j], gencov_arr[i,j], M_annot, df_prodr2.columns)
                rg_arr[j,i] = rg_arr[i,j]
                
        #save class members
        self.gencov_arr = gencov_arr
        self.rg_arr = rg_arr
        
        #write output
        self.write_output(args, M_annot)
        
        
    def load_annotations_data(self, args, df_prodr2, index_intersect):
    
        if args.annot is None and args.annot_chr is None:
            #create relevant data for a single annotation
            
            if args.M is None:
                raise ValueError('--M must be used when not using --annot or --annot-chr')
            if args.not_M_5_50 is not None:
                raise ValueError('--not-M-5-50 cannot be used without using --annot or --annot-chr')
            if args.fit_intercept:
                raise ValueError('--fit-intercept cannot be used without using --annot or --annot-chr')
            M_annot = np.ones(1) * args.M
            df_annotations_sumstats_noneg = pd.DataFrame(np.ones(len(index_intersect)), index=index_intersect, columns=['base'])
            df_sync = pd.DataFrame(index=['base'])
            df_sync['min_annot'] = 0
            df_sync['M2_5_50'] = args.M
            df_sync['is_continuous'] = False
            df_overlap = pd.DataFrame(index=['base'])
            df_overlap['base'] = args.M
            df_l2 = None
            df_w_ld = None
            
        else:
        
            #load M_annot
            if args.not_M_5_50: M_suffix = 'l2.M'
            else: M_suffix = 'l2.M_5_50'
            df_M_annot = pcgc_utils.load_dfs(args.annot, args.annot_chr, M_suffix, 'M', 'annot', header=None,use_tqdm=False)
            M_annot = df_M_annot.sum(axis=0).values
            if M_annot.shape[0] != df_prodr2.shape[1]:
                raise ValueError('.M files have a different number of columns than .prodr2 files')            
        
            #read df_sync and overlap matrix
            if args.sync is None:
                raise ValueError('--sync not provided')
            df_sync = pd.read_table(args.sync+'sync', index_col='Category')
            overlap_suffix = 'overlap' if args.not_M_5_50 else 'overlap_5_50'
            df_overlap = pd.read_table(args.sync+overlap_suffix, sep='\s+', index_col='Category')
            if df_sync.shape[0] != df_prodr2.shape[1] or not np.all(df_sync.index == df_prodr2.columns):
                raise ValueError('sync and prodr2 files must have the same annotations')
            if df_overlap.shape[0] != df_prodr2.shape[1] or not np.all(df_overlap.index == df_prodr2.columns):
                raise ValueError('overlap_matrix and prodr2 files must have the same annotations')
            
            #read SNP data files
            df_annotations_sumstats_noneg = pcgc_utils.load_dfs(args.annot, args.annot_chr, 'annot.gz', 'annot', 'annot', index_col='SNP', index_intersect=index_intersect, use_tqdm=True)
            df_annotations_sumstats_noneg.drop(columns=['CHR', 'CM', 'BP'], inplace=True)
            if df_annotations_sumstats_noneg.shape[1] != df_prodr2.shape[1] or not np.all(df_annotations_sumstats_noneg.columns == df_prodr2.columns):
                raise ValueError('annotation and prodr2 files must have the same annotations')
            df_annotations_sumstats_noneg -= df_sync['min_annot'].values
            df_list = [df_annotations_sumstats_noneg]
            if args.fit_intercept:
                df_l2 = pcgc_utils.load_dfs(args.annot, args.annot_chr, 'l2.ldscore.gz', 'l2.ldscore', 'annot', index_col='SNP', usecols=['SNP', 'baseL2'], index_intersect=index_intersect)
                df_w_ld = pcgc_utils.load_dfs(args.w_ld, args.w_ld_chr, 'l2.ldscore.gz', 'l2.ldscore', 'w-ld', index_col='SNP', usecols=['SNP', 'L2'], index_intersect=index_intersect)
                df_l2 = df_l2['baseL2']
                df_w_ld = df_w_ld['L2']
                df_list += [df_l2, df_w_ld]
            else:
                df_l2, df_w_ld = None, None
                
            #make sure that all the dfs are perfectly aligned
            for df in df_list:
                assert not df.index.duplicated().any()
                index_intersect_df = index_intersect.intersection(df.index)
                if len(index_intersect_df) < len(index_intersect):
                    raise ValueError('not all SNPs found in the annotation or LD score files - this shouldn''t happen')
                is_same = (df.index == index_intersect).all()
                if not is_same:
                    df = df.loc[index_intersect]
                    
        #restrict annotations if requested
        if args.keep_anno is not None or args.remove_anno is not None:
            #Find the set of annotations to select
            category_names = df_prodr2.columns
            remove_anno = set([] if (args.remove_anno is None) else args.remove_anno.split(','))
            keep_anno = set(category_names if (args.keep_anno is None) else args.keep_anno.split(','))
            if len(keep_anno.intersection(set(category_names))) < len(keep_anno):
                raise ValueError('-keep-anno includes non-existing annotations')
            if len(remove_anno.intersection(set(category_names))) < len(remove_anno):
                raise ValueError('-remove-anno includes non-existing annotations')
            
            #keep only the selected annotations
            anno_arr = [c for c in category_names if ((c in keep_anno) and (c not in remove_anno))]
            if len(anno_arr) < len(category_names):
                selected_anno = np.isin(category_names, anno_arr)
                M_annot = M_annot[selected_anno]
                df_prodr2 = df_prodr2.loc[anno_arr, anno_arr]
                df_annotations_sumstats_noneg = df_annotations_sumstats_noneg[anno_arr]
                df_sync = df_sync.loc[anno_arr]
                df_overlap = df_overlap.loc[anno_arr, anno_arr]
                logging.info('%d annotations remained after applying --keep-anno and --remove anno'%(df_prodr2.shape[1]))
                            
        return M_annot, df_prodr2, df_annotations_sumstats_noneg, df_sync, df_overlap, df_l2, df_w_ld

        
                
    def create_study_objects(self, args, sum_l2, category_names):
                
        #create PCGC_data objects
        if args.sumstats is None and args.sumstats_chr is None:
            raise ValueError('you muse use either --sumstats or --sumstats-chr')
        if (args.sumstats is not None) and (args.sumstats_chr is not None):
            raise ValueError('you muse use either --sumstats or --sumstats-chr')
        sumstats_prefix_list = \
            args.sumstats.split(',') if (args.sumstats is not None) else args.sumstats_chr.split(',')
        pcgc_data_list = []
        
        for f in sumstats_prefix_list:
            if args.sumstats is not None:
                pcgc_study_obj = SPCGC_Data(args, f, None, category_names, sum_l2)
            else:
                pcgc_study_obj = SPCGC_Data(args, None, f, category_names, sum_l2)
            pcgc_data_list.append(pcgc_study_obj)
            
        return pcgc_data_list
        
    
    def sync_data_files(self, args, pcgc_data_list, index_intersect):
    
        sumstats_prefix_list = \
            args.sumstats.split(',') if (args.sumstats is not None) else args.sumstats_chr.split(',')    

        #sync all the sumstats in all the files
        if index_intersect is None:
            index_intersect = pcgc_data_list[0].df_sumstats.index
        n_max = len(index_intersect)
        for pcgc_data_obj in pcgc_data_list:
            index_intersect = pcgc_data_obj.df_sumstats.index.intersection(index_intersect)
        if len(index_intersect) < MIN_NUM_SNPS:
            raise ValueError('less than %d SNPs were found in common in all the sumstats files'%(MIN_NUM_SNPS))
        if len(index_intersect) < n_max:
            logging.warning('%d SNPs are found in the annotation files and in all the sumstats files'%(len(index_intersect)))
            for pcgc_data_obj in pcgc_data_list:
                pcgc_data_obj.df_sumstats = pcgc_data_obj.df_sumstats.loc[index_intersect]
        
        #make sure that the minor alleles are synced across all files        
        if len(pcgc_data_list) > 1:
            df_alleles = pcgc_data_list[0].df_sumstats
            try:
                allele1_col = pcgc_utils.find_df_column(df_alleles, ['ALLELE2', 'A2'])
                allele0_col = pcgc_utils.find_df_column(df_alleles, ['ALLELE1', 'A1'])
            except ValueError:
                allele0_col = pcgc_utils.find_df_column(df_alleles, ['ALLELE0', 'A0'])
                allele1_col = pcgc_utils.find_df_column(df_alleles, ['ALLELE1', 'A1'])
            df_alleles = df_alleles[[allele0_col, allele1_col]]            
            for obj_i ,pcgc_data_obj in enumerate(pcgc_data_list[1:]):
                is_flipped = ((df_alleles != pcgc_data_obj.df_sumstats[[allele0_col, allele1_col]]).sum(axis=1))>0
                if np.any(is_flipped):
                    new_z = pcgc_data_obj.df_sumstats['pcgc_sumstat'].values.copy()
                    new_z[is_flipped] *= -1
                    pcgc_data_obj.df_sumstats['pcgc_sumstat'] = new_z                    
                    fname1, fname2 = sumstats_prefix_list[0], sumstats_prefix_list[obj_i]
                    logging.warning('Flipping %d SNPs in %s to match the minor alleles of %s'%(is_flipped.sum(), fname1, fname2))

        assert not index_intersect.duplicated().any()
        return index_intersect
        
        
    def write_output(self, args, M_annot):    
        
        if args.sumstats is not None:
            fname_list = args.sumstats.split(',')
        else:
            fname_list = args.sumstats_chr.split(',')    
        
        #write h2 and enrichment results to files
        for i in range(len(fname_list)):
            fname_prefix = args.out + '.' + os.path.basename(fname_list[i])
            if fname_prefix[-1] == '.': fname_prefix = fname_prefix[:-1]
            self.gencov_arr[i,i].df_enrichment.to_csv(fname_prefix+'.results', sep='\t', index=False, float_format='%0.3e', na_rep='nan')
            self.gencov_arr[i,i].df_h2.to_csv(fname_prefix+'.output', sep='\t', index=False, float_format='%0.3e')
            
            #print delete values if requested
            if args.print_delete_vals:
                np.savetxt(fname_prefix+'.delete', self.gencov_arr[i,i].delete_values.dot(M_annot))
                np.savetxt(fname_prefix+'.part_delete', self.gencov_arr[i,i].delete_values)
            
        #write rg results to file
        if len(fname_list) > 1:
            rg_arr_for_df = np.empty((len(fname_list), len(fname_list)), dtype=object)
            for i in range(len(fname_list)):
                for j in range(i+1):
                    rg_arr_for_df[i,j] = '%0.4f (%0.4f)'%(self.rg_arr[i,j].rg, self.rg_arr[i,j].rg_se)
                    rg_arr_for_df[j,i] = rg_arr_for_df[i,j]
            df_rg = pd.DataFrame(rg_arr_for_df, index=fname_list, columns=fname_list)
            rg_fname = args.out +'.rg'
            df_rg.to_csv(rg_fname, sep='\t')
                
        #write detailed rg results to file
        if args.rg_annot:
            for i in range(len(fname_list)):
                for j in range(len(fname_list)):
                    if i==j: continue
                    fname1 = os.path.basename(fname_list[i])
                    fname2 = os.path.basename(fname_list[j])
                    if fname1[-1] == '.': fname1 = fname1[:-1]
                    if fname2[-1] == '.': fname2 = fname2[:-1]
                    fname_rg = '%s.%s.%s.rg_annot'%(args.out, fname1, fname2)
                    self.rg_arr[i,j].df_annot_rg.to_csv(fname_rg, sep='\t', index=False, float_format='%0.3e')
        
                
    def create_cov_obj(self, args, oi, oj,
             df_annotations_sumstats_noneg, df_prodr2, df_sync, df_overlap, M_annot,
             df_l2, df_w_ld): 
                       
        #estimate taus and their ses
        coef, delete_values, intercept, delete_intercepts = \
            self.compute_taus(args, oi, oj,
                               df_annotations_sumstats_noneg,
                               df_prodr2,
                               df_sync,
                               M_annot,
                               df_l2, df_w_ld)
                    
        #create variance/covariance object
        M_tot = df_overlap.iloc[0,0]
        cov_obj = SPCGC_Cov(coef, delete_values, intercept, delete_intercepts, M_annot, 
                    df_prodr2.columns, df_overlap.values, M_tot, oi.var_t, oj.var_t,
                    df_sync['is_continuous'])
        return cov_obj
        
                
    def load_prodr2(self, args):
        df_prodr2 = pcgc_utils.load_dfs(args.prodr2, args.prodr2_chr, 'prodr2', 'prodr2', 'prodr2', use_tqdm=False)
        df_prodr2.index.name = 'Category'        
        df_prodr2 = df_prodr2.groupby(by=['Category']).sum()
        assert df_prodr2.shape[0] == len(df_prodr2.columns.intersection(df_prodr2.index))
        df_prodr2 = df_prodr2.loc[df_prodr2.columns, df_prodr2.columns]
        assert (df_prodr2.columns == df_prodr2.index).all()

        if args.annot is None and args.annot_chr is None:
            if df_prodr2.shape[1] > 1:
                logging.warning('Using only the first annotation in prodr2 file!!!')
                df_prodr2 = df_prodr2.iloc[:1,:1]
        
        return df_prodr2
        
    
    def load_all_snp_indices(self, args):
        index_list = []
    
        if args.annot is not None or args.annot_chr is not None:
            df_annot_index = pcgc_utils.load_dfs(args.annot, args.annot_chr, 'annot.gz', 'annot', 'annot', index_col='SNP', usecols=['SNP'], allow_duplicates=True)
            index_list.append(df_annot_index.index)
        if args.fit_intercept:
            df_l2_index = pcgc_utils.load_dfs(args.annot, args.annot_chr, 'l2.ldscore.gz', 'l2.ldscore', 'annot', index_col='SNP', usecols=['SNP'], allow_duplicates=True)
            index_list.append(df_l2_index.index)
        if args.w_ld is not None or args.w_ld_chr is not None:
            df_w_ld_index = pcgc_utils.load_dfs(args.w_ld, args.w_ld_chr, 'l2.ldscore.gz', 'l2.ldscore', 'w-ld', index_col='SNP', usecols=['SNP'], allow_duplicates=True)
            index_list.append(df_w_ld_index.index)
            
        if len(index_list) == 0:
            index_intersect = None
        elif len(index_list) == 1:
            index_intersect = index_list[0]
        else:
            index_intersect = reduce(lambda i1,i2: i1.intersection(i2), index_list)
            
        if index_intersect is not None:
            assert not index_intersect.duplicated().any()
        return index_intersect
            

    def compute_taus(self,
                   args,
                   o1, o2,
                   df_annotations_sumstats_noneg,
                   df_prodr2,
                   df_sync,
                   M_annot,
                   df_l2, df_w_ld):


        #compute the inner product of Gty1, Gty2
        df_Gty1 = o1.df_Gty
        df_Gty2 = o2.df_Gty
        intersect_columns = df_Gty1.columns.intersection(df_prodr2.columns)
        if len(intersect_columns) == 0:
            raise ValueError('Gty columns are different from the annotations names')
        if len(intersect_columns) < df_Gty1.shape[1]:
            if not np.all(df_prodr2.columns.isin(df_Gty1.columns)):
                raise ValueError('Not all annotations have Gty values')
            df_Gty1 = df_Gty1[df_prodr2.columns]
            df_Gty2 = df_Gty2[df_prodr2.columns]
        
        trace_ratios1 = o1.trace_ratios
        trace_ratios2 = o2.trace_ratios                
        if df_Gty1 is None or df_Gty2 is None or df_Gty1.shape[0] == 0 or args.fit_intercept:
            Gty12 = 0        
        else:
            #intersect df_Gty1 and df_Gty2 if required
            is_same = (df_Gty1.shape[0] == df_Gty2.shape[0]) \
                    and (df_Gty1.index == df_Gty2.index).all()
            if not is_same:
                index_intersect = df_Gty1.index.intersection(df_Gty2.index)
                df_Gty1 = df_Gty1.loc[index_intersect]
                df_Gty2 = df_Gty2.loc[index_intersect]
                logging.info('%d individuals found in both Gty files'%(df_Gty1.shape[0]))
            Gty12 = np.einsum('ij,ij->j', df_Gty1.values, df_Gty2.values) * np.sqrt(trace_ratios1*trace_ratios2)

        #fit intercept
        n1 = o1.N
        n2 = o2.N
        N = np.sqrt(n1*n2)
        sumstats1 = o1.df_sumstats['pcgc_sumstat'].values
        sumstats2 = o2.df_sumstats['pcgc_sumstat'].values
        if args.fit_intercept:
            assert df_l2 is not None
            assert df_w_ld is not None
            chi2 = sumstats1 * sumstats2 * np.sqrt(trace_ratios1 * trace_ratios2) / N
            
            #LDSC weights (uses code adapted from LDSC)
            hsq = ((chi2.mean()-1) * M_annot[0] / N) / df_l2.mean()
            hsq = np.clip(hsq, 0, 1)
            ld = np.fmax(df_l2.values, 1.0)
            w_ld = np.fmax(df_w_ld.values, 1.0)
            c = hsq * N / M_annot[0]
            intercept_temp = 1.0
            het_w = 1.0 / (2 * (intercept_temp + c*ld)**2)
            oc_w = 1.0 / w_ld
            w = np.sqrt(het_w * oc_w)
            w /= w.sum()
            
            intercept_X = np.row_stack((df_l2.values, np.ones(len(sumstats1)))).T * w[:,np.newaxis]
            intercept_Y = chi2 * w
            #intercept_Y2 = (chi2-1) * w
            intercept_XTX = intercept_X.T.dot(intercept_X)
            intercept_XTY = intercept_Y.dot(intercept_X)
            intercept = np.linalg.solve(intercept_XTX, intercept_XTY)[1] * N
        else:
            intercept = Gty12
            
            
        #compute quantities required for PCGC numerator
        z1_anno = df_annotations_sumstats_noneg.values * sumstats1[:, np.newaxis] * np.sqrt(trace_ratios1)
        z2_anno = df_annotations_sumstats_noneg.values * sumstats2[:, np.newaxis] * np.sqrt(trace_ratios2)
        z12 = np.einsum('ij,ij->j', z1_anno, z2_anno)

        #compute Z.T.dot(Y) (numer) and ZTZ (denom)        
        M_annot_sumstats2 = np.einsum('ij,ij->j', df_annotations_sumstats_noneg, df_annotations_sumstats_noneg)
        ZTZ = df_prodr2.values * \
            (np.outer(trace_ratios1,trace_ratios2)  /  np.outer(M_annot_sumstats2, M_annot_sumstats2))
        # # # M_annot_noneg2 = df_sync['M2_noneg'].values if args.not_M_5_50 else df_sync['M2_5_50_noneg'].values
        # # # ZTZ = df_prodr2.values * \
            # # # (np.outer(trace_ratios1,trace_ratios2)  /  np.outer(M_annot_sumstats2, M_annot_noneg2)); print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        ZTY = z12 / M_annot_sumstats2 - intercept
        
        #compute taus
        taus = np.linalg.solve(ZTZ, ZTY)
                
        #perform jackknife
        separators = np.floor(np.linspace(0, len(z1_anno), args.n_blocks+1)).astype(int)
        delete_values = np.empty((args.n_blocks, ZTY.shape[0]))
        delete_intercepts = np.empty(args.n_blocks)
        for block_i in range(args.n_blocks):
            b_slice = slice(separators[block_i], separators[block_i+1])
            b_M = np.einsum('ij,ij->j', df_annotations_sumstats_noneg.values[b_slice], df_annotations_sumstats_noneg.values[b_slice])
            z12_noblock_i = z12 - np.einsum('ij,ij->j', z1_anno[b_slice,:], z2_anno[b_slice,:])
            if args.fit_intercept:
                intercept_X_block = intercept_X[b_slice]
                intercept_Y_block = intercept_Y[b_slice]
                intercept_XTX_del = intercept_XTX - intercept_X_block.T.dot(intercept_X_block)
                intercept_XTY_del = intercept_XTY - intercept_Y_block.dot(intercept_X_block)
                intercept_del = np.linalg.solve(intercept_XTX_del, intercept_XTY_del)[1] * N
                delete_intercepts[block_i] = intercept_del
            else:
                intercept_del = intercept
                delete_intercepts[block_i] = intercept[0]
            ZTY_del = z12_noblock_i / (M_annot_sumstats2-b_M) - intercept_del
            delete_values[block_i] = np.linalg.solve(ZTZ, ZTY_del)
        
        #multiply by h2_factor
        deflation_ratio1 = o1.deflation_ratio
        deflation_ratio2 = o1.deflation_ratio
        deflation_ratio = np.sqrt(deflation_ratio1*deflation_ratio2)
        min_annot = df_sync['min_annot'].values
        h2_factor = 1.0 / (deflation_ratio * n1*n2 * o1.mean_Q*o2.mean_Q * (M_annot - M_annot[0]*min_annot))
        taus *= h2_factor
        delete_values *= h2_factor
        
        #scale intercept
        delete_intercepts /= N
        if args.fit_intercept: intercept /= N            
        else: intercept = intercept[0] / N
        
        #correct for annotations with negative numbers
        taus[0] -= taus[1:].dot(min_annot[1:])
        delete_values[:,0] -= np.einsum('ij,j->i', delete_values[:,1:], min_annot[1:])
        
        return taus, delete_values, intercept, delete_intercepts


    

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--sumstats', default=None, help='A comma-separated list of file prefixes (that were created by create_pcgc_sumstats.py script)')
    parser.add_argument('--sumstats-chr', default=None, help='prefix for multi-chromosome pcgc files (that were created by create_pcgc_sumstats.py script)')
    parser.add_argument('--sync', default=None, help='prefix of PCGC sync_file, created by pcgc_sync.py')
    parser.add_argument('--out', required=True, help='prefix of output file(s)')        
    parser.add_argument('--annot', default=None, help='prefix for LDSC annotation file')
    parser.add_argument('--annot-chr', default=None, help='prefix for LDSC multi-chromosome annotation file')
    parser.add_argument('--prodr2', default=None, help='prefix for r^2 product files')
    parser.add_argument('--prodr2-chr', default=None, help='prefix for multi-chromosome prefix for r^2 product files')
    parser.add_argument('--chisq-max', default=None, help='SNPs with chi^2 values above this cutoff will not be considered')
    parser.add_argument('--n-blocks', default=200, type=int, help='number of jackknife blocks')
    parser.add_argument('--not-M-5-50', default=None, action='store_true', help='If set, all reference panel SNPs will be used to estimate h2 (including ones with MAF<0.05)')
    parser.add_argument('--keep-anno', default=None, help='optional comma-separated list of annotations to use')
    parser.add_argument('--remove-anno', default=None, help='optional comma-separated list of annotations to remove')
    
    parser.add_argument('--M', default=None, type=int, help='Specify number of (common) SNPs in reference panel (not only SNPs with summary statistics). This flag can only be used when there are no annotations')
    
    parser.add_argument('--fit-intercept', default=False, action='store_true', help='fit an intercept (not recommended for PCGC)')
    parser.add_argument('--w-ld', default=None, help='LDSC weights file (only required if you want to fit an intercept)')
    parser.add_argument('--w-ld-chr', default=None, help='LDSC weights file, in multi-chromosome format (only required if you want to fit an intercept)')
    
    
    parser.add_argument('--no-Gty', default=False, action='store_true', help='Tells PCGC to assume that there are no overlapping individuals (only relevant for genetic correlation estimation)')
    parser.add_argument('--he', default=False, action='store_true', help='Use HE instead of PCGC')
    parser.add_argument('--rg-annot', default=False, action='store_true', help='This will print out a per-annotation rg table for every pair of datasets')
    
    parser.add_argument('--print-delete-vals', default=False, action='store_true', help='This will print out the jackknife delete values (may be useful for downstream scripts)')
    
    args = parser.parse_args()
    
    #print splash screen
    splash_screen()
        
    #check that the output directory exists
    if os.path.isabs(args.out) and not os.path.exists(os.path.dirname(args.out)):    
        raise ValueError('output directory %s doesn''t exist'%(os.path.dirname(args.out)))
    
    #configure logger
    pcgc_utils.configure_logger(args.out)

    #create the PCGC object (which does everything and writes output files)
    pcgc_obj = SPCGC(args)
        
