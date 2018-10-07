import numpy as np; np.set_printoptions(precision=3, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200) 
import sys
import tempfile
import subprocess
import os
import glob

def run_command(cmd):
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = proc.communicate()
    if stderr is not None:
        print '%s stderr:'%(' '.join(cmd))
        print stderr
        raise IOError()
    rc = proc.returncode
    if (rc != 0):
        print 'Command %s returned %d with the following stdout:'%(' '.join(cmd), rc)
        print stdout
        raise IOError()        
    

def test_sumstats(tmpdir):
    current_path = os.path.dirname(os.path.realpath(__file__))
    example_dir = os.path.join(current_path, 'example')
    example_files = ['example/s3','example/s2','example/s1']
    model_fname = os.path.join(example_dir, 'model')+'.'
    output_prefix = os.path.join(str(tmpdir), 'results')
    
    # #override temp dir (for testing purposes only)
    # tmpdir = '/tmp/w9Q2Im'
    # output_prefix = os.path.join(str(tmpdir), 'results')
    
    print 'temporary directory: %s'%(tmpdir)
    
    #create sumstats files
    sumstats_list = []
    study_names = []
    for example in example_files:
        fname = os.path.splitext(example)[0]
        study_name = os.path.splitext(os.path.basename(fname))[0]
        study_names.append(study_name)
        print 'Creating summary statistics for study %s...'%(study_name)
        sumstats_list.append(os.path.join(str(tmpdir), '%s.'%(study_name)))
        for chr_num in xrange(1,23):
            model_chr_fname = os.path.join(example_dir, 'model')+'.%d.'%(chr_num)
            out_fname = os.path.join(str(tmpdir), '%s.%d'%(study_name, chr_num))
            #out_fname = os.path.join(example_dir, 'results', '%s.%d'%(study_name, chr_num))
            ss_command = ['python', 'pcgc_sumstats_creator.py']
            ss_command += ['--bfile', fname]
            ss_command += ['--covar', fname+'.cov']
            ss_command += ['--prev', '0.01']
            ss_command += ['--pheno', fname+'.phe']
            ss_command += ['--frqfile', model_chr_fname]
            ss_command += ['--annot', model_chr_fname]
            ss_command += ['--out', out_fname]
            ss_command += ['--sync', model_fname]
            run_command(ss_command)
            
    #run PCGC
    pcgc_command = ['python', 'pcgc_main.py']
    pcgc_command += ['--frqfile-chr', model_fname]
    pcgc_command += ['--annot-chr', model_fname]
    pcgc_command += ['--prodr2-chr', model_fname]
    pcgc_command += ['--sumstats-chr', ','.join(sumstats_list)]
    pcgc_command += ['--out', output_prefix]
    pcgc_command += ['--sync', model_fname]
    run_command(pcgc_command)
    
    #check output files
    dir_gold = os.path.join(example_dir, 'results')
    dir_test = str(tmpdir)
    
    #check rg
    df_rg_gold = pd.read_table(os.path.join(dir_gold, 'results.rg'), index_col=0)
    df_rg_test = pd.read_table(os.path.join(dir_test, 'results.rg'), index_col=0)
    #remove the parentheses
    df_rg_gold.columns = np.arange(len(df_rg_gold.columns))
    df_rg_test.columns = np.arange(len(df_rg_gold.columns))
    df_rg_gold.index = np.arange(len(df_rg_gold.index))
    df_rg_test.index = np.arange(len(df_rg_gold.index))
    for c in df_rg_gold.columns:
        df_rg_gold[c] = df_rg_gold[c].str.replace(r"\(.*\)","")
        df_rg_test[c] = df_rg_test[c].str.replace(r"\(.*\)","")
    df_rg_gold = df_rg_gold.astype(np.float)
    df_rg_test = df_rg_test.astype(np.float)
    assert np.allclose(df_rg_gold, df_rg_test)
    
    #check main output fields
    for study_name in study_names:
        df_output_test = pd.read_table(os.path.join(dir_test, 'results.%s.output'%(study_name)), index_col='Quantity')
        df_output_gold = pd.read_table(os.path.join(dir_gold, 'results.%s.output'%(study_name)), index_col='Quantity')
        assert np.allclose(df_output_test, df_output_gold)

    #check enrichment
    for study_name in study_names:
        df_enr_test = pd.read_table(os.path.join(dir_test, 'results.%s.results'%(study_name)), index_col='Category')
        df_enr_gold = pd.read_table(os.path.join(dir_gold, 'results.%s.results'%(study_name)), index_col='Category')
        assert np.allclose(df_enr_test, df_enr_gold, equal_nan=True)
    
                

if __name__ == '__main__':
    temp_dir = os.path.join(tempfile._get_default_tempdir(), next(tempfile._get_candidate_names()))
    os.mkdir(temp_dir)
    test_sumstats(temp_dir)