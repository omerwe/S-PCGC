import numpy as np; np.set_printoptions(precision=3, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200) 
import scipy.stats as stats
import scipy.linalg as la
import sys
import time
import tempfile
import itertools
import subprocess
import os
from tqdm import tqdm
import string
import random

PLINK_EXE = '/home/ow14/plink/plink'

class CC_Study:

    def __init__(self, mafs, h2, beta, beta_covariates, prev, n, df_map, Z_shared_con,
                 frac_cases=0.5, use_liab=False, plink_exe=None):
    
        #generate SNPs and liabilities for an entire population
        m = len(mafs)
        num_gen = np.maximum(int(float(n) / float(prev)), 25000)        
        Z_snps =  (np.random.random(size=(num_gen, m)) < mafs).astype(np.int)
        Z_snps += (np.random.random(size=(num_gen, m)) < mafs).astype(np.int)        
        if Z_shared_con is not None:
            Z_snps = np.concatenate((Z_shared_con, Z_snps), axis=0)
            num_gen = Z_snps.shape[0]
        
        snp_std = np.sqrt(2*mafs*(1-mafs))
        liab = Z_snps.dot(beta/snp_std) - 2*mafs.dot(beta/snp_std)
        liab += np.random.randn(num_gen) * np.sqrt(1-h2)    #add noise        
        
        #generate covariates and update liabilities
        if beta_covariates is None or len(beta_covariates)==0:
            X = np.empty((num_gen,0))
        else:
            X = np.random.randint(2, size=(num_gen, len(beta_covariates)))
            liab += X.dot(beta_covariates)
            
        #find cases and controls
        affection_cutoff = np.percentile(liab, 100*(1-prev))
        cases = (liab>affection_cutoff)
        controls = ~cases

        #define the phenotype y
        if use_liab: y = liab            
        else: y = cases.astype(np.int)
        
        #generate random string for study name
        random_str = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6))
        
        #select cases
        num_cases = int(n*frac_cases)
        assert cases.sum() >= n/2.0, 'not enough cases found'
        assert controls.sum() >= n/2.0, 'not enough controls found'
        case_inds = np.where(liab >= affection_cutoff)[0]
        assert len(case_inds >= num_cases)
        case_inds = np.random.permutation(case_inds)[:num_cases]
        case_iid = np.array(['case%d_%s'%(i+1,random_str) for i in xrange(len(case_inds))])
        
        #select controls        
        num_controls = int(int(n*(1-frac_cases)))
        if Z_shared_con is None:
            control_inds1 = np.zeros(0, dtype=np.int)
            num_shared = 0
        else:
            num_shared = Z_shared_con.shape[0]
            control_inds1 = np.arange(num_shared)
        control1_iid = np.array(['shared_control%d'%(i+1) for i in xrange(len(control_inds1))])
        if len(control_inds1) >= num_controls:
            randperm = np.random.permutation(len(control_inds1))[:num_controls]
            control_inds = control_inds1[randperm]
            control_iid = control1_iid[randperm]
        else:
            num_controls2 = num_controls - len(control_inds1)
            is_control2 = liab < affection_cutoff
            is_control2[:num_shared] = False
            control_inds2 = np.where(is_control2)[0]
            assert len(control_inds2 >= num_controls2)
            control_inds2 = np.random.permutation(control_inds2)[:num_controls2]
            control2_iid = np.array(['control%d_%s'%(i+1,random_str) for i in xrange(len(control_inds2))])
            control_inds = np.concatenate((control_inds1, control_inds2))
            control_iid = np.concatenate((control1_iid, control2_iid))
        
        #save selected sample
        study_inds = np.concatenate((case_inds, control_inds))
        study_iid = np.concatenate((case_iid, control_iid))
        self.Z_snps = Z_snps[study_inds]
        self.y = y[study_inds]
        self.X = X[study_inds]

        #normalize SNPs
        snps_mean = 2*mafs
        snps_std = np.sqrt(2*mafs*(1-mafs))
        #snps_mean = Z_snps.mean(axis=0)
        #maf_est = snps_mean/2.0
        #snps_std = Z_snps.std(axis=0)
        #snps_std = np.sqrt(2*maf_est*(1-maf_est))        
        self.Z = (self.Z_snps - snps_mean) / snps_std
        self.Z_r = self.Z
        self.X_PCs = np.empty((self.Z_snps.shape[0], 0))
        self.X_all = self.X
        self.study_iid = study_iid
        
        #save plink executable
        if plink_exe is None: self.plink_exe = PLINK_EXE
        else: self.plink_exe = plink_exe
        
        #save df_map and prev
        self.df_map = df_map
        self.prev = prev
        
        
    def compute_PCs(self, num_PCs):
    
        if num_PCs==0:
            return
        
        K = self.Z.dot(self.Z.T) / self.Z.shape[1]
        s,U = la.eigh(K)        
        X_PCs = U[:,-num_PCs:]
    
        #regress PCs out of Z            
        mean_diag = np.diag(K).mean()
        Z_r = self.Z - X_PCs.dot(np.linalg.solve(X_PCs.T.dot(X_PCs), X_PCs.T.dot(self.Z)))
        
        #apply trace correction (to ensure that the mean trace value is 1.0)
        n, m = self.Z_snps.shape[0], self.Z_snps.shape[1]
        mean_diag_r = np.einsum('ij,ij->', Z_r, Z_r) / (n*m)
        Z_r *= np.sqrt(mean_diag / mean_diag_r)
        
        self.Z_r = Z_r
        self.X_PCs = X_PCs
        self.eigvals = s
        
        self.X_all = np.concatenate((self.X, self.X_PCs), axis=1)
        
        
    def write_plink_file(self):

        #create a temporary file name
        plink_fname = os.path.join(tempfile._get_default_tempdir(), next(tempfile._get_candidate_names()))
        iid = self.study_iid
        
        #create a phenotypes file
        df_pheno = pd.DataFrame(self.y, columns=['pheno'])
        df_pheno.insert(0, 'fid', iid)
        df_pheno.insert(1, 'iid', iid)
        df_pheno.to_csv(plink_fname+'.phe', sep='\t', index=False, header=True)
        
        #create a covariates file if required        
        if self.X_all.shape[1] > 0:            
            cov_names = ['cov%d'%(i+1) for i in xrange(self.X_all.shape[1])]
            df_cov = pd.DataFrame(self.X_all, columns=cov_names)
            df_cov.insert(0, 'fid', iid)
            df_cov.insert(1, 'iid', iid)
            df_cov.to_csv(plink_fname+'.cov', sep='\t', index=False, header=True)         
        
        #create a (text) plink map file                
        self.df_map.to_csv(plink_fname+'.map', sep='\t', index=False, header=False)
        
        #Create a df of phased SNPs
        assert list(np.unique(self.Z_snps)) == [0,1,2]
        n, m = self.Z_snps.shape[0], self.Z_snps.shape[1]
        Z_haploid = np.empty((self.Z_snps.shape[0], 2*m), dtype=np.int)
        is_mat_first = np.random.random(size=(n,m)) < 0.5
        Z_p = np.zeros((n,m))
        Z_m = np.zeros((n,m))
        Z_p[self.Z_snps==2] = 1
        Z_m[self.Z_snps==2] = 1
        Z_m[(self.Z_snps==1) & is_mat_first] = 1
        Z_p[(self.Z_snps==1) & (~is_mat_first)] = 1
        Z_haploid[:, ::2] = Z_p
        Z_haploid[:, 1::2] = Z_m
        snp_names = ['snp%d'%(i+1) for i in xrange(m)]
        snp_col_names = list(itertools.chain.from_iterable([(s+'_1', s+'_2') for s in snp_names]))
        df_Z_haploid = pd.DataFrame(Z_haploid+1, columns=snp_col_names)
        
        #create a (text) plink ped file
        df_ped = pd.DataFrame(iid, columns=['fid'])
        df_ped['iid'] = iid
        df_ped['father'] = 0
        df_ped['mother'] = 0
        df_ped['sex'] = 1
        df_ped['pheno'] = self.y
        df_ped = pd.concat((df_ped, df_Z_haploid), axis=1)
        df_ped.to_csv(plink_fname+'.ped', sep='\t', index=False, header=False)
            
        #create a binary plink file
        cmdLine = [self.plink_exe, '--file', plink_fname, '--make-bed', '--out', plink_fname]
        proc = subprocess.Popen(cmdLine, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = proc.communicate()
        if (stderr is not None):
            print 'plink stderr:'
            print stderr
            raise Exception()        
        assert os.path.exists(plink_fname+'.bed')
            
        #create pve files
        if self.X_PCs.shape[1]>0:
            pve = self.eigvals / self.eigvals.sum()
            df_pve = pd.DataFrame(pve, columns=['pve'])
            df_pve.to_csv(plink_fname+'.pve', header=None, index=None)
        
        self.plink_fname = plink_fname

            
        
class CC_Studies:
    def __init__(self, annotations, h2_arr, beta_covar_arr, m, n_arr, c_arr, prev_arr, h2_c_arr, 
                 num_shared_con, frac_cases_arr=None, use_liab_arr=None, plink_exe=None):
    
        #default setting is 50% cases
        num_studies = len(n_arr)
        if frac_cases_arr is None:
           frac_cases_arr = np.array([0.5] * num_studies)
        if use_liab_arr is None:
            use_liab_arr = np.zeros(num_studies, dtype=np.bool)

        #create chromosome numbers
        chr_arr = np.ones(m, dtype=np.int)
        bin_size = (m / 2) / 21
        bin_i = m / 2
        for chr_num in xrange(2,23):
            chr_arr[bin_i : bin_i+bin_size] = chr_num
            bin_i += bin_size
        if bin_i < len(chr_arr)-1:
            chr_arr[bin_i:]=22
        self.chr_arr = chr_arr
            
        #create a df of plink map file
        self.create_plink_maf(m)            
        
        #save class members
        self.annotations = annotations
        self.h2_arr = h2_arr
        self.beta_covar_arr = beta_covar_arr
        self.m = m
        self.n_arr = n_arr
        self.c_arr = c_arr
        self.prev_arr = prev_arr
        self.h2_c_arr = h2_c_arr
        self.num_shared_con = num_shared_con
        self.frac_cases_arr = frac_cases_arr
        self.use_liab_arr = use_liab_arr
        self.plink_exe = plink_exe
        
        
    def simulate_studies(self):
        
        #generate MAFs
        mafs = 0.05 + np.random.random(size=self.m) * 0.45
        
        #generate betas (beta_covar_arr has dimensions [anno_i, study_i, study_j])
        num_studies = len(self.n_arr)
        num_anno = self.annotations.shape[1]
        beta3d = np.empty((self.m, num_studies, num_anno))
        for anno_i in xrange(num_anno):
            assert len(np.unique(np.sign(np.diag(self.beta_covar_arr[anno_i]))))==1
            L_sign = np.sign(self.beta_covar_arr[anno_i,0,0])
            L = la.cholesky(self.beta_covar_arr[anno_i] * L_sign, lower=True)            
            assert np.allclose(L.dot(L.T)*L_sign, self.beta_covar_arr[anno_i])
            #print 'generating betas for annotation %d'%(anno_i+1)
            for m_i in tqdm(range(self.m), disable=True):
                beta3d[m_i, :, anno_i] = L.dot(np.random.randn(num_studies))*L_sign
            
        #Compute betas of each SNP ii each study according to its annotations
        beta_arr = np.empty((num_studies, self.m))
        for study_i in xrange(num_studies):
            beta_arr[study_i] = np.einsum('ij,ij->i', beta3d[:,study_i,:], self.annotations)
                    
        #generate beta_covariates
        beta_covariates = []
        for study_i in xrange(num_studies):
            if self.c_arr[study_i]==0:
                beta_covariates_i = np.array([])
            else:
                beta_covariates_i = np.random.randn(self.c_arr[study_i]) * np.sqrt(self.h2_c_arr[study_i] / self.c_arr[study_i])
            beta_covariates.append(beta_covariates_i)        
        
        #generate shared controls
        if self.num_shared_con==0:
            Z_shared_con = None
        else:
            Z_shared_con =  (np.random.random(size=(self.num_shared_con, self.m)) < mafs).astype(np.float)
            Z_shared_con += (np.random.random(size=(self.num_shared_con, self.m)) < mafs).astype(np.float)                    
        
        #create studies
        studies_arr = []
        for i in xrange(num_studies):
            print 'creating study %d/%d'%(i+1, num_studies)
            study = CC_Study(mafs, self.h2_arr[i],
                    beta_arr[i], beta_covariates[i], self.prev_arr[i], self.n_arr[i],                    
                    df_map=self.df_map, Z_shared_con=Z_shared_con, use_liab = self.use_liab_arr[i],
                    frac_cases=self.frac_cases_arr[i], plink_exe=self.plink_exe)
            studies_arr.append(study)
            
        #save class members        
        self.studies_arr = studies_arr
        self.mafs = mafs

    def write_files(self, multi_chrom):
        #write ref files
        print 'writing annotation files...'
        self.write_ref_files(multi_chrom)
        
        #write plink files
        print 'writing plink files...'
        for s in self.studies_arr:
            s.write_plink_file()
                 
        
    
    def get_plink_fnames(self):
        plink_fnames = [s.plink_fname for s in self.studies_arr]
        return plink_fnames
        
        
    def create_plink_maf(self, m):
        snp_names = ['snp%d'%(i+1) for i in xrange(m)]        
        df_map = pd.DataFrame(snp_names, columns=['SNP'])            
        df_map.insert(0, 'CHR', self.chr_arr)
        df_map['CM'] = np.arange(m) / 10.0
        df_map['BP'] = np.arange(m) * 10        
        self.df_map = df_map
        
        
    def write_ref_files(self, multi_chrom):
    
        ref_fname = os.path.join(tempfile._get_default_tempdir(), next(tempfile._get_candidate_names()))
    
        #create MAF files
        df_mafs = self.df_map[['CHR', 'SNP']].copy()        
        df_mafs['MAF'] = self.mafs
        df_mafs['A0'] = 2
        df_mafs['A1'] = 1
        df_mafs.to_csv(ref_fname+'.frq', sep='\t', index=False, header=True)
        for chr_num in self.chr_arr:
            anno_chr = self.annotations[self.chr_arr==chr_num]
            df_mafs_chr = df_mafs.query('CHR == %d'%(chr_num))
            df_mafs_chr.to_csv(ref_fname+'.%d.frq'%(chr_num), sep='\t', index=False, header=True)
            
        #create a sync file
        annotation_names = ['anno_%d'%(anno_i+1) for anno_i in xrange(self.annotations.shape[1])]        
        min_annot = np.min(self.annotations, axis=0)
        min_annot[min_annot>0]=0
        df_sync = pd.Series(min_annot, index=annotation_names)
        df_sync.index.name = 'Category'
        df_sync.to_csv(ref_fname+'.sync', sep='\t', float_format='%0.5e')        

        #create prod_r^2 files and M_annot files
        df_prod_r2 = pd.DataFrame(((self.annotations-min_annot)**2).T.dot((self.annotations-min_annot)**2), index=annotation_names, columns=annotation_names)
        df_M = pd.DataFrame(np.row_stack(self.annotations.sum(axis=0)).T, columns=annotation_names)
        if not multi_chrom:
            df_prod_r2.to_csv(ref_fname+'.prodr2', sep='\t', index=True, header=True)
            df_M.to_csv(ref_fname+'.l2.M_5_50', header=False, index=False, float_format='%0.3f', sep='\t')
        else:
            for chr_num in self.chr_arr:
                anno_chr = self.annotations[self.chr_arr==chr_num]
                df_prod_r2_chr = pd.DataFrame((anno_chr**2).T.dot(anno_chr**2), index=annotation_names, columns=annotation_names)
                df_prod_r2_chr.to_csv(ref_fname+'.%d.prodr2'%(chr_num), sep='\t', index=True, header=True)
                df_M_chr = pd.DataFrame(np.row_stack(anno_chr.sum(axis=0)).T, columns=annotation_names)
                df_M_chr.to_csv(ref_fname+'.%d.l2.M_5_50'%(chr_num), header=False, index=False, float_format='%0.3f', sep='\t')
        
        #create annotation files
        df_anno = self.df_map[['CHR', 'BP', 'SNP', 'CM']].copy()
        for anno_i, anno_name in enumerate(annotation_names):
            df_anno[anno_name] = self.annotations[:, anno_i]
        if not multi_chrom:
            df_anno.to_csv(ref_fname+'.annot.gz', sep='\t', index=False, header=True, compression='gzip')
        else:
            for chr_num in self.chr_arr:
                df_anno_chr = df_anno.query('CHR == %d'%(chr_num))
                df_anno_chr.to_csv(ref_fname+'.%d.annot.gz'%(chr_num), sep='\t', index=False, header=True, compression='gzip')        
            
        #save files prefix
        self.ref_fname = ref_fname
        
        

    
        

        
        
        

