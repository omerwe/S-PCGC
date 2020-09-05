import numpy as np; np.set_printoptions(precision=3, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200) 
import scipy.stats as stats
import scipy.linalg as la
import sys
import time
from sklearn.linear_model import LogisticRegression
import tempfile
import itertools
import subprocess
import os
import logging
import imp
logging.getLogger().setLevel(logging.CRITICAL)

import pcgc_sumstats_creator
import pcgc_main
import pcgc_simulator
import pcgc_r2

PLINK_EXE = '/home/ow14/plink/plink'
GCTA_EXE = '/home/ow14/apps/gcta_1.26.0/gcta64'

class Dummy:
    def __init__(self):
        pass
        
        
def run_gcta(gcta_exe, plink_fname, prev):
    ''' Estimate h2 using gcta GREML
        gcta_exe: executable for gcta
        plink_exe: executable for plink        
    '''
    
    #create a GRM
    cmdLine = [gcta_exe, '--bfile', plink_fname, '--make-grm-bin', '--out', plink_fname]    
    proc = subprocess.Popen(cmdLine, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = proc.communicate()
    if (stderr is not None):
        print('GCTA make-grm error:')
        print(stderr)
        raise Exception()
    assert os.path.exists(plink_fname+'.grm.bin')
    
    #run GREML analysis
    cmdLine = [gcta_exe, '--grm-bin', plink_fname, '--reml', '--pheno', plink_fname+'.phe', '--prevalence', str(prev), '--out', plink_fname]
    if os.path.exists(plink_fname + '.cov'):
        cmdLine += ['--qcovar', plink_fname+'.cov']
    proc = subprocess.Popen(cmdLine, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = proc.communicate()
    if (stderr is not None):
        print('GCTA REML error:')
        print(stderr)
        raise Exception()
    rc = proc.returncode
    if (rc != 0):
        print('Command %s returned %d with the following stdout:'%(' '.join(cmdLine), rc))
        print(stdout)
        raise IOError()        
    assert os.path.exists(plink_fname+'.hsq')
    
    #read gcta results
    df_reml = pd.read_table(plink_fname+'.hsq', sep='\s+', skiprows=[5,6], index_col='Source')
    gcta_h2_est = df_reml.loc['V(G)/Vp_L', 'Variance']
    
    return gcta_h2_est
    
        

def compute_Q(prev, X=None, y=None):    
    
    P = np.mean(y>y.mean())
    if (X is None) or X.shape[1]==0:
        tau = stats.norm(0,1).isf(prev)
        phi_tau = stats.norm(0,1).pdf(tau)
        Q = P*(1-P) / (prev**2 * (1-prev)**2) * phi_tau**2
        return Q
    else:
        assert y is not None
        
    logreg = LogisticRegression(penalty='l2', C=500000, fit_intercept=True)         
    logreg.fit(X, y)
    Pi = logreg.predict_proba(X)[:,1]
        
    K = prev
    Ki = K*(1-P) / (P*(1-K)) * Pi / (1 + K*(1-P) / (P*(1-K))*Pi - Pi)
    tau_i = stats.norm(0,1).isf(Ki)
    tau_i[Ki>=1.] = -999999999
    tau_i[Ki<=0.] = 999999999   
    phi_tau_i = stats.norm(0,1).pdf(tau_i)
    
    u_prefix = phi_tau_i / np.sqrt(Pi*(1-Pi)) / (Ki + (1-Ki)*(K*(1-P))/(P*(1-K)))
    u0 = u_prefix * K*(1-P) / (P*(1-K)) * Pi
    u1 = u_prefix * (1-Pi)
    Q = np.outer(u0,u0) + np.outer(u0,u1) + np.outer(u1,u0) + np.outer(u1,u1)    
    assert np.all(~np.isnan(Q))
    
    return Q


def pcgc(K_list, y1, y2, prev1, prev2, X1, X2, remove_diag=True):
    
    #compute Q and standardize y
    if len(np.unique(y1)) == 2:
        Q1 = compute_Q(prev1, X1, y2)
        P1 = np.mean(y1>y1.mean())
        y1_norm = (y1-P1) / np.sqrt(P1*(1-P1))
    else:
        Q1 = 1
        y1_norm = (y1-y1.mean()) / y1.std()
        
    if len(np.unique(y2)) == 2:
        Q2 = compute_Q(prev2, X2, y2)
        P2 = np.mean(y2>y2.mean())
        y2_norm = (y2-P2) / np.sqrt(P2*(1-P2))
    else:
        Q2 = 1
        y2_norm = (y2-y2.mean()) / y2.std()
        
    #compute shared Q
    Q = np.sqrt(Q1*Q2)

    #create V (equivalent to denom)
    num_anno = len(K_list)
    V = np.empty((num_anno,num_anno))
    for i in range(num_anno):
        KQ_i = K_list[i]*Q
        for j in range(i,num_anno):
            KQ_j = K_list[j]*Q            
            V[i,j] = np.einsum('ij,ij->', KQ_i, KQ_j)
            if remove_diag:
                V[i,j] -= np.diag(KQ_i).dot(np.diag(KQ_j))
            V[j,i] = V[i,j]
            
    #create R (equivalent to numer)
    R = np.empty(num_anno)
    yy = np.outer(y1_norm,y2_norm)
    for i in range(num_anno):
        KQ_i = K_list[i]*Q
        R[i] = np.einsum('ij,ij->', yy, KQ_i)
        if remove_diag:
            R[i] -= np.diag(KQ_i).dot(y1_norm*y2_norm)
        
    pcgc_est = np.linalg.solve(V, R)    
    return pcgc_est
    
    

    
def create_chr_extract_file(plink_fname, chr_num):
    extract_fname = os.path.join(tempfile._get_default_tempdir(), next(tempfile._get_candidate_names()))
    df_bim = pd.read_table(plink_fname+'.bim', sep='\s+', usecols=[0,1])
    df_bim.columns = ['CHR', 'SNP']
    df_bim = df_bim.query('CHR==%d'%(chr_num))
    assert df_bim.shape[0] > 0
    df_bim[['SNP']].to_csv(extract_fname, header=False, index=False)
    return extract_fname
        
def run_plink_linreg(plink_exe, plink_fname, out_fname, n, chr_num=None):

    #create plink command
    plink_command = [plink_exe, '--bfile', plink_fname, '--allow-no-sex', '--assoc', '--out', out_fname, '--pheno', plink_fname+'.phe']
    if os.path.exists(plink_fname+'.cov'):
        plink_command += ['--covar', plink_fname+'.cov']
    if chr_num is not None:
        extract_fname = create_chr_extract_file(plink_fname, chr_num)
        plink_command += ['--extract', extract_fname]

    #run plink command
    proc = subprocess.Popen(plink_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = proc.communicate()
    if (stderr is not None):
        print('plink stderr:')
        print(stderr)
        raise Exception()        
    assert os.path.exists(out_fname+'.qassoc')
    
    #add some required fields and save to the correct file name
    df_linreg = pd.read_table(out_fname+'.qassoc', sep='\s+')
    df_linreg['N'] = n
    df_linreg['Z'] = stats.norm(0,1).isf(df_linreg['P'] / 2.0) * np.sign(df_linreg['BETA'])
    df_linreg['A0'] = 2
    df_linreg['A1'] = 1
    df_linreg.to_csv(out_fname, sep='\t', index=False)
    
    
def run_pcgc_sumstats_creator(prev, ref_fname, plink_fname, out_fname, study_obj, multi_chrom, chr_num=None):
    #create an args object
    sumstats_args = Dummy()
    sumstats_args.prev = prev
    sumstats_args.bfile = plink_fname
    sumstats_args.bgen = None
    sumstats_args.pheno = plink_fname+'.phe'
    if os.path.exists(plink_fname+'.cov'):
        sumstats_args.covar = plink_fname+'.cov'
    else:
        sumstats_args.covar = None
    sumstats_args.pheno_col = 1
    sumstats_args.keep = None
    sumstats_args.remove = None
    sumstats_args.exclude = None
    sumstats_args.frqfile = ref_fname+'.'
    sumstats_args.frqfile_chr = None
    sumstats_args.sync = ref_fname+'.'
    sumstats_args.min_maf = None
    
    if multi_chrom:
        sumstats_args.annot = None
        sumstats_args.annot_chr = ref_fname+'.'
    else:
        sumstats_args.annot = ref_fname+'.'
        sumstats_args.annot_chr = None
    if os.path.exists(plink_fname+'.pve'):
        sumstats_args.pve = plink_fname+'.pve'
        num_c = study_obj.X.shape[1]
        num_PCs = study_obj.X_PCs.shape[1]
        sumstats_args.covars_regress = ','.join(['cov%d'%(i) for i in range(num_c+1, num_PCs+1)])
    else:
        sumstats_args.covars_regress = None
        sumstats_args.pve = None
        
    if chr_num is None:
        sumstats_args.extract = None
    else:
        sumstats_args.extract = create_chr_extract_file(plink_fname, chr_num)    
        
    #create a sumstats object and write its output to a temp file
    imp.reload(pcgc_sumstats_creator); sumstats_creator = pcgc_sumstats_creator.PCGC_Sumstats(args=sumstats_args)
    sumstats_creator.compute_all_sumstats()
    sumstats_creator.write_output(out_fname)        
        
     
def create_sumstats(studies_obj, multi_chrom):
    
    sumstats_prefixes = []    
    ref_fname = studies_obj.ref_fname    
    for study_i in range(num_studies):
    
        #extract the data of study i
        study_obj = studies_obj.studies_arr[study_i]
        plink_fname = study_obj.plink_fname
        prev = study_obj.prev

        #create a file name for the temporary summary statistics
        ss_fname = os.path.join(tempfile._get_default_tempdir(), next(tempfile._get_candidate_names()))
        sumstats_prefixes.append(ss_fname)        
        
        #if not a case-control study, run Plink
        if len(np.unique(study_obj.y)) > 2:
            n = len(study_obj.y)
            if multi_chrom:
                for chr_num in range(1,23):
                    run_plink_linreg(PLINK_EXE, plink_fname, ss_fname+'.%d'%(chr_num), n, chr_num=chr_num)
            else:
                run_plink_linreg(PLINK_EXE, plink_fname, ss_fname, n, chr_num=None)
                
        #if it's a case-control study
        else:
            if multi_chrom:
                for chr_num in range(1,23):
                    run_pcgc_sumstats_creator(prev, ref_fname, plink_fname, ss_fname+'.%d'%(chr_num), study_obj, multi_chrom, chr_num=chr_num)
            else:
                run_pcgc_sumstats_creator(prev, ref_fname, plink_fname, ss_fname, study_obj, multi_chrom, chr_num=None)
        
    return sumstats_prefixes
    
    
    
def run_pcgc_sumstats(ref_fname, sumstats_prefixes, multi_chrom, use_he):

    out_fname_prefix = os.path.join(tempfile._get_default_tempdir(), next(tempfile._get_candidate_names()))
        
    #create an args object
    spcgc_args = Dummy()    
    if multi_chrom:
        spcgc_args.sumstats_chr = ','.join([c+'.' for c in sumstats_prefixes])
        spcgc_args.sumstats = None
        spcgc_args.annot_chr = ref_fname+'.'
        spcgc_args.annot = None
        spcgc_args.prodr2_chr = ref_fname+'.'
        spcgc_args.prodr2 = None
        spcgc_args.frqfile_chr = ref_fname+'.'
        spcgc_args.frqfile = None
    else:
        spcgc_args.sumstats = ','.join([c+'.' for c in sumstats_prefixes])
        spcgc_args.sumstats_chr = None
        spcgc_args.annot = ref_fname+'.'
        spcgc_args.annot_chr = None
        spcgc_args.prodr2 = ref_fname+'.'
        spcgc_args.prodr2_chr = None
        spcgc_args.frqfile = ref_fname+'.'
        spcgc_args.frqfile_chr = None
    spcgc_args.out = out_fname_prefix
    spcgc_args.no_Gty = False
    spcgc_args.not_M_5_50 = False
    spcgc_args.n_blocks = 200
    spcgc_args.chisq_max = None
    spcgc_args.he = use_he
    spcgc_args.fit_intercept = False
    spcgc_args.no_annot = False
    spcgc_args.rg_annot = True
    spcgc_args.keep_anno = None
    spcgc_args.print_delete_vals = True
    spcgc_args.remove_anno = None
    spcgc_args.sync = ref_fname+'.'
    
    #run S-PCGC
    imp.reload(pcgc_main); pcgc_obj = pcgc_main.SPCGC(args=spcgc_args)
    return pcgc_obj, out_fname_prefix
    
def create_effects_covariance_matrix(annotations, h2_arr, rg_matrix):

    #create the marginal variances of each annotation in each study
    num_anno = annotations.shape[1]
    num_studies = len(h2_arr)
    sigma2_anno_base = np.abs(np.random.randn(num_anno))
    #sigma2_anno_base[1]=1e-6; print 'sigma2_anno_base[1]=1e-6!!!!!!!!!!!!'
    snp_var_base = annotations.dot(sigma2_anno_base)        
    sigma2_anno_arr = np.empty((num_studies, num_anno))
    for study_i in range(num_studies):
        sigma2_anno_arr[study_i] = sigma2_anno_base * h2_arr[study_i] / snp_var_base.sum()    
        if not np.all(annotations.dot(sigma2_anno_arr[study_i]) > 0):
            raise ValueError()
        
    #create the covariance matrix
    beta_covar_arr = np.empty((num_anno, num_studies, num_studies))    
    for anno_i in range(num_anno):
        beta_covar_arr[anno_i] = np.diag(sigma2_anno_arr[:, anno_i])
        b = beta_covar_arr[anno_i]
        for s_i in range(num_studies):
            for s_j in range(s_i):                
                b[s_i,s_j] = rg_matrix[s_i,s_j] * np.sqrt(b[s_i,s_i] * b[s_j,s_j])
                b[s_j,s_i] = b[s_i,s_j]
    assert np.all(~np.isnan(beta_covar_arr)) 

    return beta_covar_arr
    
    
    
    
    
def pcgc_direct_all(studies_obj, annotations, use_PCs):
    num_studies = len(studies_obj.studies_arr)
    num_anno = annotations.shape[1]
    
    #Prepare quantities for correction of annotations with negative numbers
    M = annotations.sum(axis=0).astype(np.float)
    min_annot = np.min(annotations, axis=0); min_annot[min_annot>0]=0
    
    #create all GRMs
    K_arr = np.empty((num_anno, num_studies, num_studies), dtype=np.object)
    for anno_i in range(num_anno):
        w = annotations[:,anno_i] - min_annot[anno_i]; assert np.all(w >= 0)
        for study_i in range(num_studies):
            o_i = studies_obj.studies_arr[study_i] 
            Z_i = (o_i.Z_r if use_PCs else o_i.Z)
            for study_j in range(study_i, num_studies):
                o_j = studies_obj.studies_arr[study_j]
                Z_j = (o_j.Z_r if use_PCs else o_j.Z)
                K = (Z_i * w**2).dot(Z_j.T) / w.dot(w)
                
                #zero out entries of identical individuals
                if study_i != study_j:
                    df_iid_i = pd.DataFrame(np.arange(len(o_i.study_iid)), index=o_i.study_iid, columns=['ind_pos'])
                    df_iid_j = pd.DataFrame(np.arange(len(o_j.study_iid)), index=o_j.study_iid, columns=['ind_pos'])
                    df_iid_both = df_iid_i.merge(df_iid_j, left_index=True, right_index=True, suffixes=('_i', '_j'))
                    print('Studies %d,%d include %d shared individuals'%(study_i+1, study_j+1, df_iid_both.shape[0]))
                    for ind, row in df_iid_both.iterrows():
                        ind_pos_i, ind_pos_j = row['ind_pos_i'], row['ind_pos_j']
                        K[ind_pos_i,ind_pos_j]=0
                       
                #save the GRM
                K_arr[anno_i, study_i, study_j] = K
                K_arr[anno_i, study_j, study_i] = K
                
    #compute h2 / genetic covariance
    gencov = np.zeros((num_anno, num_studies, num_studies))
    for study_i in range(num_studies):
        y_i = studies_obj.studies_arr[study_i].y
        X_i = studies_obj.studies_arr[study_i].X_all
        prev_i = prev_arr[study_i]
        for study_j in range(study_i+1):
            y_j = studies_obj.studies_arr[study_j].y
            X_j = studies_obj.studies_arr[study_j].X_all
            prev_j = prev_arr[study_j]
            pcgc_est = pcgc(K_arr[:,study_i,study_j], y_i, y_j, prev_i, prev_j, X_i, X_j, remove_diag=(study_i==study_j))
            
            #apply corretion for annotations with negative numbers
            tau = pcgc_est / (M-annotations.shape[0]*min_annot)
            tau[0] -= tau[1:].dot(min_annot[1:])
            pcgc_est = tau*M
            
            gencov[:, study_i, study_j] = pcgc_est
            gencov[:, study_j, study_i] = pcgc_est
    
    #print out all the results
    for study_i in range(num_studies):
        print('h2 study %d: %0.3f'%(study_i+1, gencov[:, study_i, study_i].sum()), gencov[:, study_i, study_i])
    for study_i in range(num_studies):
        for study_j in range(study_i+1, num_studies):
            print('gencor studies %d,%d:'%(study_i+1, study_j+1), end=' ')
            print('%0.3f'%(gencov[:, study_i, study_j].sum() / np.sqrt(gencov[:, study_i, study_i].sum()*gencov[:, study_j, study_j].sum())))
    
    
    
    
    
def simulate_studies(num_snps, h2_arr, prev_arr, n_arr, anno_freq, rg=0.5, rg_matrix=None,
                    c_arr=None, h2_c_arr=None, num_shared_con=0, use_liab=False,
                    num_pcs_arr=None):
    '''
        Simulate a set of studies with functional annotations
    
        num_snps: number of SNPs
        h2_arr: An array of true simulated h2 for each study
        prev_arr: An array of disease prevalance for each study
        n_arr: An array of sample sizes for each study
        anno_freq: The %SNPs in each annotation. The first number must be 1.0
        rg: The overall genetic correlation between each pair of studies
        rg_matrix: A matrix with the rg between each pair of studies. If this is provided, rg will be ignored
        c_arr: An array with the number of binary covariates in each study
        h2_c_arr: An array with the liability variance explained by the covariates in each study
        num_shared_con: The number of controls shared between these studies
        use_liab: If True, the liability is observed (i.e., the phenotype is quantitative)
        num_pcs_arr: The number of PCs to use in each study
    '''
        
    num_studies = len(h2_arr)
    m = num_snps
    num_anno = len(anno_freq)
    
    #Replace Nones with default values
    if c_arr is None:
        c_arr = np.zeros(num_studies, dtype=np.int)        
    if h2_c_arr is None:
        h2_c_arr = np.ones(num_studies) * 0.25
    if num_pcs_arr is None:
        num_pcs_arr = np.zeros(num_studies, dtype=np.int)        
    
    #generate annotations
    anno_freq = [1,0.5]#,0.25,0.25]     #the number of SNPs in each annotations
    assert anno_freq[0]==1    
    annotations = (np.random.random(size=(m,num_anno))<=anno_freq).astype(np.float)
    #annotations[:,1] = np.ones(m); annotations[m/2:, 1] = -1
    #annotations[:,1] = np.random.randn(m)
    assert np.all(annotations[:,0] == 1)
    
    #generate an rg matrix is not provided
    if rg_matrix is None:
        rg_matrix = np.ones((num_studies, num_studies))*rg
    
    #generate covariance matrices of every annotations in every pair of studies
    for try_num in range(100):
        try:
            beta_covar_arr = create_effects_covariance_matrix(annotations, h2_arr, rg_matrix)
        except ValueError:
            continue
        break
    
    #create simulated datasets
    imp.reload(pcgc_simulator); studies_obj = pcgc_simulator.CC_Studies(annotations=annotations, h2_arr=h2_arr, beta_covar_arr=beta_covar_arr, m=m, n_arr=n_arr, c_arr=c_arr, prev_arr=prev_arr, h2_c_arr=h2_c_arr, num_shared_con=num_shared_con, use_liab_arr=np.ones(num_studies)*int(use_liab)); studies_obj.simulate_studies()

    #apply PC correction
    for study_i in range(num_studies):
        studies_obj.studies_arr[study_i].compute_PCs(num_pcs_arr[study_i])
        
    #run direct PCGC-s
    # pcgc_direct_all(studies_obj, annotations, use_PCs=False)
    # print
    # pcgc_direct_all(studies_obj, annotations, use_PCs=True)
    # print
    
    #create plink and annotation files
    studies_obj.write_files(multi_chrom=False)
    
    return studies_obj
    
    
def pcgc_h2(studies_obj, use_he=False):
    sumstats_prefixes = create_sumstats(studies_obj, multi_chrom=False)
    pcgc_obj, results_prefix = run_pcgc_sumstats(studies_obj.ref_fname, sumstats_prefixes, multi_chrom=False, use_he=use_he)
    h2_estimates = [pcgc_obj.gencov_arr[i,i].tot for i in range(len(studies_obj.studies_arr))]
    return h2_estimates
    
def gcta_h2(studies_obj):
    h2_gcta_list = []
    for study in studies_obj.studies_arr:
        h2_gcta = run_gcta(GCTA_EXE, study.plink_fname, study.prev)
        h2_gcta_list.append(h2_gcta)
    return np.array(h2_gcta_list)
    
    
    

if __name__ == '__main__':    
    
    h2_arr=[0.5, 0.4]                   #h2 of each study
    num_studies = len(h2_arr)           #number of studies
    prev_arr = [0.01]*num_studies       #prevalence of each study
    n_arr = [500]*num_studies          #sample sizes
    num_snps = 600                     #number of SNPs
    c_arr = [0]*num_studies             #number of covariates in each study
    h2_c_arr = [0.25]*num_studies       #liability variance explained by covariates in each study
    num_shared_con = 0                  #number of shared controls across studies
    use_liab = False                    #specifies if the liability is observed
    anno_freq = [1, 0.5]                #frequency of each annotation
    num_pcs_arr = [0]*num_studies
    
    num_experiments = 100
    gcta_bias_list = []
    pcgc_bias_list = []
    for exp_num in range(num_experiments):
        studies_obj = simulate_studies(num_snps, h2_arr, prev_arr, n_arr, anno_freq, 
                            c_arr=c_arr, h2_c_arr=h2_c_arr, num_shared_con=num_shared_con, use_liab=use_liab,
                            num_pcs_arr=num_pcs_arr)
        h2_experiment_pcgc = pcgc_h2(studies_obj, use_he=use_liab)
        h2_experiment_gcta = gcta_h2(studies_obj)
        for h2_i in range(len(h2_experiment_pcgc)):
            gcta_bias = h2_experiment_gcta[h2_i] - h2_arr[h2_i]
            gcta_bias_list.append(gcta_bias)
            pcgc_bias = h2_experiment_pcgc[h2_i] - h2_arr[h2_i]
            pcgc_bias_list.append(pcgc_bias)
            
        gcta_bias_avg = np.mean(gcta_bias_list)
        pcgc_bias_avg = np.mean(pcgc_bias_list)
        gcta_bias_std = np.std(gcta_bias_list)
        pcgc_bias_std = np.std(pcgc_bias_list)
        print('Experiment %d  PCGC bias: %0.3f (%0.3f)  GCTA bias: %0.3f (%0.3f)'%(exp_num+1, pcgc_bias_avg, pcgc_bias_std, gcta_bias_avg, gcta_bias_std))
            
        
    
    
    
    
