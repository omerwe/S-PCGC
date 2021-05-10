import numpy as np; np.set_printoptions(precision=3, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200)
import scipy.stats as stats
import os
import time
import sys
import scipy.linalg as la
from tqdm import tqdm
import logging
from pandas_plink import read_plink
#from bgen_reader import read_bgen, convert_to_dosage
from sklearn.linear_model import LogisticRegression
import pcgc_utils

MAX_SNPS_IN_MEMORY = int(5e7)

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

def splash_screen():
    print('*********************************************************************')
    print('* S-PCGC sumstats creator')
    print('* Version 2.0.0')
    print('* (C) 2018 Omer Weissbrod')
    print('*********************************************************************')
    print()


class PCGC_Sumstats:
    def __init__(self, args):
    
        #read phenotypes file
        df_pheno = self.read_pheno_file(args)
            
        #read covariates
        if args.covar is not None:
            df_covar = pd.read_table(args.covar, sep='\s+', dtype=str)
            for c in df_covar.columns[2:]:
                df_covar[c] = df_covar[c].astype(np.float)
                if np.any(df_covar[c].isnull()):
                    raise ValueError('covariate %s includes missing values. Plesae impute them or remove individuals with missing covariates'%(c))
            df_covar = self.add_fid_iid_index(df_covar)
            
            #merge individuals across phenotypes and covariates
            index_intersect = df_pheno.index.intersection(df_covar.index)
            if len(index_intersect) < len(df_pheno.index):
                if len(index_intersect)==0:
                    raise ValueError('no individuals have both both phenotype and covariates data')
                df_pheno = df_pheno.loc[index_intersect]
                df_covar = df_covar.loc[index_intersect]
                logging.info('%d individuals have both phenotypes and covariates data'%(df_covar.shape[0]))
                
                
        #read mafs file if it exists
        if args.frqfile is not None or args.frqfile_chr is not None:        
            df_maf = pcgc_utils.load_dfs(args.frqfile, args.frqfile_chr, 'frq', 'frq', 'frqfile', index_col='SNP')            
        else:
            df_maf = None
            logging.warning('MAF file not provided! We will use the in-sample MAF estimates (this is highly unrecommended!)')
                
        #read bfile or bgen file
        if args.bfile is not None:
            assert args.bgen is None, '--bfile and --bgen cannot both be specified'            
            self.bfile, df_pheno, df_maf, self.num_snps, self.sample_size = self.read_plink(args, df_pheno, df_maf)
            self.genetic_format = 'plink'
        elif args.bgen is not None:
            assert args.bfile is None, '--bfile and --bgen cannot both be specified'
            raise NotImplementedError('bgen functionality not yet implemented')
            self.genetic_format = 'bgen'
        else:
            raise ValueError('either --bfile or --bgen must be specified')
            
        #save MAFs
        if df_maf is not None:
            maf_col = self.find_df_column(df_maf, ['MAF', 'FRQ', 'FREQ', 'A1Freq'], 'MAF')
            self.mafs = df_maf[maf_col]
        else:
            self.mafs = None
            
        #Extract relevant covariates and compute the Cholesky factorization of the hat matrix    
        if args.covar is None:
            C = None
            self.C_regress = None
            self.L_CTC = None
            covars_regress_cols = []
        else:
            #reorder individuals to make sure that the covariates match the other files
            is_same = (df_covar.shape[0] == df_pheno.shape[0]) and (df_covar.index == df_pheno.index).all()
            if not is_same:
                df_covar = df_covar.loc[df_pheno.index]
            C = df_covar.iloc[:, 2:]
            
            #extract relevant covariates
            if args.covars_regress is None:
                self.C_regress = None
                self.L_CTC = None
                covars_regress_cols = []
            else:
                assert args.covar is not None
                covars_regress_cols = args.covars_regress.split(',')
                for c in covars_regress_cols:
                    if c not in df_covar.columns:
                        raise ValueError('%s is not in the covariates file'%(c))
                self.C_regress = df_covar[covars_regress_cols].values
                if np.any(np.isnan(self.C_regress)):
                    raise ValueError('Covariates include missing values')
                
                if np.linalg.matrix_rank(self.C_regress) < self.C_regress.shape[1]:
                    raise ValueError('Some of the covariates are linearly dependent on other covariates --- please double check your covariates input')

                #compute the Cholesky factorization of the hat matrix
                self.L_CTC = la.cho_factor(self.C_regress.T.dot(self.C_regress))
                if np.any(np.isnan(self.L_CTC[1])):
                    raise ValueError('Covariates matrix could not be factorized --- please check if your covariates are linearly dependent')
                
        
        #load pve file if it exists
        if args.pve is None:
            if args.covars_regress is not None:
                raise ValueError('--covars_regress must be used with --pve')
            self.pve = []
        else:
            df_pve = pd.read_table(args.pve, header=None)            
            if df_pve.shape[1] > 1:
                raise ValueError('pve file must include only a single column')
            if df_pve.shape[0] < len(covars_regress_cols):
                raise ValueError('There are fewer pve values than covariates to regress')
            df_pve.columns = ['pve']
            if df_pve.shape[0] > len(covars_regress_cols):
                logging.warning('There are more pve values than covariates to regress. Using only top %d pve values'%(len(covars_regress_cols)))
                df_pve = df_pve.sort_values('pve', ascending=False).head(len(covars_regress_cols))
            self.pve = df_pve.iloc[:,0].values
            
                            
        #transform phenotypes to 0/1
        y = df_pheno.iloc[:,-1]
        num_pheno = len(np.unique(y))
        if num_pheno == 1:
            raise ValueError('only one phenotype value found!')
        elif num_pheno > 2:
            raise ValueError('phenotype file must include only cases and controls')
        y = (y>y.mean()).astype(np.int).values
        
        #compute PCGC statistics
        P = y.mean()        
        u0, u1, ty, self.var_t = self.compute_pcgc_stats(args.prev, P, y, C)
        y_norm = (y-P) / np.sqrt(P*(1-P))
        self.z_coeff = y_norm * (u0+u1)
        self.mean_Q = np.mean((u0 + u1)**2)
        self.prev = args.prev
        
        #load annotations
        self.load_annotations(args.annot, args.annot_chr, args.sync)

        
    def load_annotations(self, anno, anno_chr, sync_prefix):
    
        #extract SNP names
        if self.genetic_format == 'plink':
            index_snpnames = self.bfile['df_bim'].index
        elif self.genetic_format == 'bgen':
            raise NotImplementedError()
        else:
            raise ValueError('illegal option')            
            
        #load annotations (or create baseline annotation only)
        if (anno is None) and (anno_chr is None):
            df_annotations = pd.DataFrame(np.ones(len(index_snpnames)), index=index_snpnames, columns=['base'])            
        else:
            df_annotations = pcgc_utils.load_dfs(anno, anno_chr, 'annot.gz', 'annot', 'annot', index_col='SNP')
            df_annotations.drop(columns=['CHR', 'CM', 'BP'], inplace=True)
            if not np.allclose(df_annotations.iloc[:,0], 1):
                raise ValueError('The first annotation must be the base annotation (all SNPs must have value 1.0)')
                
            #apply min_annot correction to ensure no negative values
            category_names = df_annotations.columns
            if sync_prefix is None:
                raise ValueError('--annot and --annot-chr must be used together with --sync')
            df_sync = pd.read_table(sync_prefix+'sync', index_col='Category')
            if df_sync.shape[0] != len(category_names) or not np.all(df_sync.index == category_names):
                raise ValueError('Annotations in sync file do not match those in annotations/prodr2 files')
            min_annot = df_sync['min_annot'].values
            df_annotations -= min_annot
                
        #remove annotations for unused SNPs
        is_same = (df_annotations.shape[0] == len(index_snpnames)) and \
                  (df_annotations.index == index_snpnames).all()
        if not is_same:
            has_anno = (index_snpnames.isin(df_annotations.index)).all()
            if not has_anno:
                raise ValueError('not all SNPs have annotations')
            df_annotations = df_annotations.loc[index_snpnames]
        
        #save the df in a class member
        self.df_annotations_noneg = df_annotations
            
            

        
        
    def write_output(self, out_prefix):
        logging.info('Writing output files...')
    
        #write sumstats file
        if self.genetic_format == 'plink':
            df_sumstats = self.bfile['df_bim'][['snp', 'a0', 'a1']].copy()
            df_sumstats.columns = ['SNP', 'A0', 'A1']
            df_fam = self.bfile['df_fam']
        elif self.genetic_format == 'bgen':
            raise NotImplementedError()
        else:
            raise ValueError('illegal option')        
        df_sumstats['pcgc_sumstat'] = self.z
        df_sumstats['N'] = self.sumstats_N
        df_sumstats.to_csv(out_prefix + '.sumstats.gz', index=False, sep='\t', float_format='%0.6e', compression='gzip')
        
        #write Gty file
        Gty = np.sqrt(self.G_diag) * self.z_coeff[:, np.newaxis]
        df_Gty = pd.DataFrame(Gty, index=df_fam.index, columns=self.df_annotations_noneg.columns)
        df_Gty.insert(0, 'iid', df_fam.iloc[:,1])
        df_Gty.insert(0, 'fid', df_fam.iloc[:,0])        
        df_Gty.to_csv(out_prefix + '.Gty.gz', index=False, sep='\t', float_format='%0.6e', compression='gzip')
        
        #write diagGRM files
        df_diagGRM = df_fam.iloc[:,:2].copy()
        for c_i, c in enumerate(self.df_annotations_noneg.columns):
            df_diagGRM['diag_G_deflate_%s'%(c)] = self.G_diag[:,c_i]
            df_diagGRM['diag_G_nodeflate_%s'%(c)] = self.G_noregress_diag[:,c_i]
        df_diagGRM.to_csv(out_prefix + '.diagGRM.gz', index=False, sep='\t', float_format='%0.6e', compression='gzip')

        #write other stats file        
        df_dicts = []
        df_dicts.append({'Property':'N', 'Value':self.sample_size})
        df_dicts.append({'Property':'prevalence', 'Value':self.prev})
        df_dicts.append({'Property':'var_t', 'Value':self.var_t})
        df_dicts.append({'Property':'mean_Q', 'Value':self.mean_Q})
        for pve_i, pve in enumerate(self.pve):
            df_dicts.append({'Property':'pve%d'%(pve_i+1), 'Value':pve})
        m_anno2 = np.einsum('ij,ij->j', self.df_annotations_noneg, self.df_annotations_noneg)
        for anno_i, anno_name in enumerate(self.df_annotations_noneg):
            df_dicts.append({'Property':'M2_%s'%(anno_name), 'Value':m_anno2[anno_i]})
        df_otherstats = pd.DataFrame(df_dicts)
        df_otherstats.to_csv(out_prefix + '.otherstats', index=False, sep='\t', float_format='%0.6e')
        
        
    def compute_all_sumstats(self, chunk_size=None):
    
        if chunk_size is None:
            chunk_size = int(float(MAX_SNPS_IN_MEMORY) / float(self.sample_size))
            logging.info('Setting chunk size to %d to keep at most %d alleles in memory at once'%(chunk_size, chunk_size*self.sample_size))
        
        self.z = np.zeros(self.num_snps)
        self.sumstats_N = np.zeros(self.num_snps)
        self.G_diag = np.zeros((self.sample_size, self.df_annotations_noneg.shape[1]))
        if self.L_CTC is not None:
            self.G_noregress_diag = np.zeros((self.sample_size, self.df_annotations_noneg.shape[1]))
        logging.info('Starting summary statistics computation...')
        for snp1 in tqdm(range(0, self.num_snps, chunk_size)):            
            snp2 = snp1 + chunk_size            
            self.set_locus(snp1, snp2)
            annot_noneg = self.df_annotations_noneg.iloc[snp1:snp2].values
            self.sumstats_N[snp1:snp2] = self.n
            if self.L_CTC is not None:
                self.G_noregress_diag += np.einsum('ij,ij,jk,jk->ik', self.X, self.X, annot_noneg, annot_noneg)
            
            #regress covariates out of the genotypes if required
            if self.L_CTC is None:
                X = self.X
            else:                
                X = self.X - self.C_regress.dot(la.cho_solve(self.L_CTC, self.C_regress.T.dot(self.X)))
            
            #compute required summary statistics
            self.z[snp1:snp2] = self.z_coeff.dot(X)            
            ###self.G_diag += (X**2).dot(annot_noneg**2)
            ###assert np.allclose((X**2).dot(self.df_annotations_noneg.iloc[snp1:snp2].values**2), np.einsum('ij,ij,jk,jk->ik', X, X, annot_noneg, annot_noneg))
            self.G_diag += np.einsum('ij,ij,jk,jk->ik', X, X, annot_noneg, annot_noneg)

        if self.L_CTC is None:
            self.G_noregress_diag = self.G_diag
            
        
    def set_locus(self, snp1, snp2):
        
        #read SNPs from bfile/bgen file
        if self.genetic_format == 'plink':            
            X = self.bfile['bed'][:, snp1:snp2].compute().astype(np.float)
        elif self.genetic_format == 'bgen':
            raise NotImplementedError('bgen functionality not yet implemented')
        else:
            raise ValueError('internal error - unknonwn genetic format')
            
        #normalize SNPs (the loop is faster than vector computations, believe it or not)
        n = X.shape[0]
        self.n = np.zeros(X.shape[1], dtype=np.int)
        for j in range(0, X.shape[1]):
            newsnp = X[:, j]
            is_miss = np.isnan(newsnp)
            num_miss = is_miss.sum()
            self.n[j] = n - num_miss
            newsnp[is_miss] = 0
            if self.mafs is None:                
                snp_avg = newsnp.sum() / float(n - num_miss)
                snp_maf = 1 - snp_avg/2.0
                ###snp_var = (newsnp.dot(newsnp)) / float(n - num_miss) - snp_avg**2
                ###assert np.all(snp_var>0)
            else:
                snp_maf = self.mafs[snp1+j]
                snp_avg = 2-2*snp_maf
            newsnp -= snp_avg
            newsnp[is_miss]=0
            newsnp /= (-np.sqrt(2*snp_maf*(1-snp_maf)))

        #save the resulting genotypes
        self.X = X
        
        
        
    def compute_pcgc_stats(self, prev, P, y, covars=None):
        
        if (covars is None):
            Pi = P
        else:
            logreg = LogisticRegression(penalty='l2', C=500000, fit_intercept=True)            
            logreg.fit(covars, y)
            Pi = logreg.predict_proba(covars)[:,1]            
            
        K = prev
        Ki = K*(1-P) / (P*(1-K)) * Pi / (1 + K*(1-P) / (P*(1-K))*Pi - Pi)
        tau_i = stats.norm(0,1).isf(Ki)
        if covars is not None:
            tau_i[Ki>=1.] = -999999999
            tau_i[Ki<=0.] = 999999999   
        phi_tau_i = stats.norm(0,1).pdf(tau_i)
        
        u_prefix = phi_tau_i / np.sqrt(Pi*(1-Pi)) / (Ki + (1-Ki)*(K*(1-P))/(P*(1-K)))
        u0 = u_prefix * K*(1-P) / (P*(1-K)) * Pi
        u1 = u_prefix * (1-Pi)
        ty = (y-Pi) / np.sqrt(Pi * (1-Pi))
        
        #compute liability variance due to covariates
        if covars is None: var_t = 0
        else: var_t = self.varLiab_covar(prev, tau_i, y)
        
        return u0, u1, ty, var_t
        
        
    #compute liability variance due to covariates
    def varLiab_covar(self, prev, tau_i, y):
        cases = y>y.mean()
        controls = ~cases
        var_E_t_given_y = prev * (1-prev) * (tau_i[cases].mean() - tau_i[controls].mean())**2
        E_var_t_given_y = prev * np.var(tau_i[cases]) + (1-prev) * np.var(tau_i[controls])
        var_t = var_E_t_given_y + E_var_t_given_y
        return var_t        
        
        
    def find_df_column(self, df, strings_to_find, df_name):
        
        if isinstance(strings_to_find, str):
            strings_to_find = [strings_to_find]
            
        is_relevant_col = np.zeros(df.shape[1], dtype=np.bool)
        for str_to_find in strings_to_find:
            is_relevant_col = is_relevant_col | (df.columns.str.upper() == str_to_find.upper())
        if np.sum(is_relevant_col)==0:
            raise ValueError('No matching column found among: %s in %s file'%(str(strings_to_find), df_name))
        elif np.sum(is_relevant_col)>1:
            raise ValueError('Too many matching columns found among: %s in %s file'%(str(strings_to_find), df_name))
        else:
            return df.columns[is_relevant_col][0]
            
            
        
    def read_pheno_file(self, args):        
        #read phenotypes from file
        df_pheno = pd.read_table(args.pheno, sep='\s+', usecols=[0,1,args.pheno_col+1], dtype=str)
        df_pheno.iloc[:,-1] = df_pheno.iloc[:,-1].astype(np.int)
        df_pheno = self.add_fid_iid_index(df_pheno)

        #apply --keep
        if args.keep is not None:
            df_keep = pd.read_table(args.keep, sep='\s+', header=None, usecols=[0,1])            
            df_keep = self.add_fid_iid_index(df_keep)
            df_pheno = df_pheno.merge(df_keep.iloc[:,:0], how='inner', left_index=True, right_index=True)
            if df_pheno.shape[0]==0:
                raise ValueError('no individuals remained after applying --keep!')
            logging.info('%d phenotyped individuals remained after applying --keep'%(df_pheno.shape[0]))
         
        #apply --remove
        if args.remove is not None:
            df_remove = pd.read_table(args.remove, sep='\s+', header=None, usecols=[0,1])
            df_remove.columns = ['fid', 'iid']
            df_remove = self.add_fid_iid_index(df_remove)            
            df_pheno = df_pheno.loc[~df_pheno.index.isin(df_remove.index)]
            if df_pheno.shape[0]==0:
                raise ValueError('no individuals remained after applying --remove!')
            logging.info('%d phenotyped individuals remained after applying --remove'%(df_pheno.shape[0]))
            
        
        #remove individuals with non-numeric phenotypes (e.g. NANs)
        is_bad_pheno = ~(df_pheno.iloc[:,-1].apply(np.isreal))
        is_bad_pheno = is_bad_pheno | df_pheno.iloc[:,-1].isnull()
        if np.any(is_bad_pheno):
            logging.warning('Removing %d individuals with non-numeric phenotypes'%(is_bad_pheno.sum()))
            df_pheno = df_pheno.loc[~is_bad_pheno]
            df_pheno.iloc[:,-1] = df_pheno.iloc[:,-1].astype(np.float)
        
        #verify that we have only cases and controls
        num_pheno = len(np.unique(df_pheno.iloc[:,-1]))
        if num_pheno == 1:
            raise ValueError('only one phenotype value found!')
        elif num_pheno > 2:
            raise ValueError('phenotype file must include only cases and controls')
        
        return df_pheno
        
        
        

    def read_plink(self, args, df_pheno, df_maf=None):
        #read plink file from disk
        logging.info('Reading plink file from disk (this may take a while...)')
        (df_bim, df_fam, bed) = read_plink(args.bfile)
        df_bim.set_index('snp', drop=False, inplace=True)
        df_bim['a0'] = df_bim['a0'].astype('str')
        df_bim['a1'] = df_bim['a1'].astype('str')
        bed = bed.T
        
        #apply --extract
        if args.extract is not None:
            df_extract = pd.read_table(args.extract, sep='\s+', header=None, usecols=[0])
            is_good_snp = df_bim['snp'].isin(df_extract.iloc[:,0])
            if not np.all(is_good_snp):
                df_bim = df_bim.loc[is_good_snp]
                bed = bed[:, is_good_snp.values]
            if df_bim.shape[0]==0:
                raise ValueError('no SNPs remained after applying --extract!')
            logging.info('%d SNPs remained after applying --extract'%(df_bim.shape[0]))
        
        #apply --exclude
        if args.exclude is not None:
            df_exclude = pd.read_table(args.exclude, sep='\s+', header=None, usecols=[0])
            is_good_snp = ~(df_bim['snp'].isin(df_exclude.iloc[:,0]))
            if not np.all(is_good_snp):
                df_bim = df_bim.loc[is_good_snp]
                bed = bed[:, is_good_snp.values]
            if df_bim.shape[0]==0:
                raise ValueError('no SNPs remained after applying --exclude!')
            logging.info('%d SNPs remained after applying --exclude'%(df_bim.shape[0]))
            
        #filter out SNPs without maf, or with a small MAF
        if df_maf is not None:
            is_good_snp = df_bim['snp'].isin(df_maf.index)
            if not np.all(is_good_snp):
                df_bim = df_bim.loc[is_good_snp]
                bed = bed[:, is_good_snp.values]
                if df_bim.shape[0] < 10:
                    raise ValueError('<10 SNPs remained after intersecting with frqfile')
                logging.info('%d SNPs remained after removing SNPs without MAFs'%(df_bim.shape[0]))
                
            #make sure that MAFs and SNPs are synchronized
            is_same = (df_maf.shape[0] == df_bim.shape[0]) and (df_maf.index == df_bim['snp']).all()
            if not is_same:
                df_maf = df_maf.loc[df_bim['snp']]                
                
            #find SNPs with different alleles between plink file and MAF file        
            try:
                allele1_col = self.find_df_column(df_maf, ['ALLELE2', 'A2'], 'MAF')
                allele0_col = self.find_df_column(df_maf, ['ALLELE1', 'A1'], 'MAF')
            except ValueError:
                allele0_col = self.find_df_column(df_maf, ['ALLELE0', 'A0'], 'MAF')
                allele1_col = self.find_df_column(df_maf, ['ALLELE1', 'A1'], 'MAF')
            df_maf[allele0_col] = df_maf[allele0_col].astype('str')
            df_maf[allele1_col] = df_maf[allele1_col].astype('str')
            is_consistent_snp = (df_maf[allele1_col] == df_bim['a1']) & (df_maf[allele0_col] == df_bim['a0'])
            is_consistent_snp = is_consistent_snp | (df_maf[allele1_col] == df_bim['a0']) & (df_maf[allele0_col] == df_bim['a1'])
            if not np.all(is_consistent_snp):
                logging.info('%d SNPs will be removed because they have different alleles in plink and MAF files'%(np.sum(~is_consistent_snp)))
                df_maf = df_maf.loc[is_consistent_snp].copy()
                df_bim = df_bim.loc[is_consistent_snp].copy()
                bed = bed[:, is_consistent_snp.values]

            #flip MAFs of flipped appeles
            is_flipped = ((df_maf[allele1_col] == df_bim['a0']) & (df_maf[allele0_col] == df_bim['a1'])).values
            maf_col = self.find_df_column(df_maf, ['MAF', 'FRQ', 'FREQ', 'A1Freq'], 'MAF')
            if np.any(is_flipped):
                new_mafs = df_maf[maf_col].values
                new_mafs[is_flipped] = 1 - new_mafs[is_flipped]
                df_maf.loc[:,maf_col] = new_mafs
                logging.warning('The minor allele of %d SNPs in the MAF file is different from the one in the plink file'%(np.sum(is_flipped)))
                
            #perform a MAF-based filtering
            if args.min_maf > 0:                
                is_good_snp = df_maf[maf_col] > args.min_maf
                if not np.all(is_good_snp):
                    df_bim = df_bim.loc[is_good_snp]
                    bed = bed[:, is_good_snp.values]
                    logging.info('%d SNPs remained after removing SNPs with MAF<%0.3f'%(df_bim.shape[0], args.min_maf))
                    
            #make sure that MAFs and SNPs are synchronized (again)
            is_same = (df_maf.shape[0] == df_bim.shape[0]) and (df_maf.index == df_bim['snp']).all()
            if not is_same:
                df_maf = df_maf.loc[df_bim['snp']]
                
        #filter individuals without a phenotype (or covariates), while preserving plink file order
        df_fam = self.add_fid_iid_index(df_fam)
        has_pheno = df_fam.index.isin(df_pheno.index)
        if not np.all(has_pheno):
            bed = bed[has_pheno, :]
            df_fam = df_fam.loc[has_pheno]
            df_pheno = df_pheno.loc[df_fam.index]
            if df_pheno.shape[0]==0:
                raise ValueError('no individuals found in both the plink and phenotype files!')            
            logging.info('%d individuals have both genotypes and phenotypes'%(df_pheno.shape[0]))
        assert np.all(df_fam.index.isin(df_pheno.index))
        df_pheno = df_pheno.loc[df_fam.index]
                            
        #check that everything is consistent
        if df_bim.shape[0] < 10:
            raise ValueError('<10 SNPs remained after all filtering stages')
        assert (df_pheno.index == df_fam.index).all()
        if df_maf is not None:
            assert (df_maf.index == df_bim.index).all()
        logging.info('%d individuals remained after all filtering stages'%(df_pheno.shape[0]))
        logging.info('%d SNPs remained after all filtering stages'%(df_bim.shape[0]))

        #create a bfile dictionary
        bfile = {'df_bim':df_bim, 'df_fam':df_fam, 'bed':bed}
        num_snps = df_bim.shape[0]
        sample_size = df_fam.shape[0]
        return bfile, df_pheno, df_maf, num_snps, sample_size


    def add_fid_iid_index(self, df):
        df['fid_iid'] = df.iloc[:,0].map(str) + '_____' + df.iloc[:,1].map(str)
        df.set_index('fid_iid', drop=True, inplace=True)
        return df
        

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--bfile', required=True, type=str, help='plink file')
    parser.add_argument('--sync', default=None, help='prefix of PCGC sync_file, created by pcgc_sync.py')
    parser.add_argument('--prev', required=True, type=float, help='trait prevalence')
    parser.add_argument('--out', required=True, type=str, help='prefix of output files')    
    parser.add_argument('--bgen', default=None, type=str, help='bgen file (not yet supported)')
    parser.add_argument('--pheno', default=None, type=str, help='phenotypes file in plink format (with header)')
    parser.add_argument('--pheno-col', default=1, type=int, help='phenotype number to use, if phenotype files includes multiple phenotype columns')
    parser.add_argument('--frqfile', default=None, help='MAFs file, in plink format')
    parser.add_argument('--frqfile-chr', default=None, help='Multi-chromosome MAFs file prefix, in plink format')
    parser.add_argument('--covar', default=None, type=str, help='covariates file, in plink format (with header)')
    parser.add_argument('--covars-regress', default=None, type=str, help='names of covariates to regress out of the genotypes (usually these will be PCs, as computed by e.g. FlashPCA2)')
    parser.add_argument('--pve', default=None, type=str, help='pve file with %%variance explained by top PCs, without a header, as created by e.g. FlashPCA2 (no header, one number per line)')
    
    parser.add_argument('--annot', help='annotations file prefix in LDSC format')
    parser.add_argument('--annot-chr', help='annotations file prefix in LDSC format, splitted across chromosomes')    
    
    parser.add_argument('--remove', default=None, type=str, help='individuals to remove')
    parser.add_argument('--keep', default=None, type=str, help='individuals to keep')
    parser.add_argument('--exclude', default=None, type=str, help='SNPs to remove')
    parser.add_argument('--extract', default=None, type=str, help='SNPs to keep')
    
    parser.add_argument('--min-maf', type=float, default=0, help='Filter out SNPs with MAF smaller than this')
    parser.add_argument('--chunk-size', type=int, default=None, help='#SNPs to load to memory at once. Larger values will enable faster analysis but will take up more memory')
    
    args = parser.parse_args()
    
    #print splash screen
    splash_screen()
    
    #check that the output directory exists
    if os.path.isabs(args.out) and not os.path.exists(os.path.dirname(args.out)):
        raise ValueError('output directory %s doesn''t exist'%(os.path.dirname(args.out)))
        
    #configure logger
    configure_logger(args.out)

    #run analysis    
    sumstats_creator = PCGC_Sumstats(args)
    sumstats_creator.compute_all_sumstats(args.chunk_size)
    sumstats_creator.write_output(args.out)
