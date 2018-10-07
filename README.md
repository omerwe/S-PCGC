# S-PCGC
Heritability, genetic correlation and functional enrichment estimation for case-control studies.

S-PCGC is an adaptation of stratified LD score regression [(S-LDSC)](https://www.nature.com/articles/ng.3404) for case-control studies. Like S-LDSC, S-PCGC can estimate genetic heritability, genetic correlation and functional enrichment. However, S-PCGC is explicitly designed for case-control studies rather than quantitative phenotypes. This can make a large differnce, especially in the presence of strong non-genetic risk factors. Parts of the S-PCGC code are adapted from S-LDSC with permission.

The main features of S-PCGC are:
1. **S-PCGC is designed for case-control studies**. Such studies include many subtleties not accounted for by methods designed for quantitative phenotypes like S-LDSC.
2. **Seamless integration with S-LDSC format**. S-PCGC accepts the same input and creates the same output files as S-LDSC.
3. **Computational efficiency**. S-PCGC can analyze datasets with hundreds of thousands of individuals and dozens of functional annotations in a few hours.

S-PCGC is described in the following paper:
<br>
[Estimating SNP-Based Heritability and Genetic Correlation in Case-Control Studies Directly and with Summary Statistics. The American Journal of Human Genetics, 103(1) 89-99.](https://www.sciencedirect.com/science/article/pii/S0002929718301952). This paper however, doesn't include a description of functional enrichment estimation. We will post a bioRxiv preprint with details soon.

<br><br>
# Installation
S-PCGC is designed for Python 2.7, and depends on the following freely available Python packages:
* [numpy](http://www.numpy.org/) and [scipy](http://www.scipy.org/)
* [scikit-learn](http://scikit-learn.org/stable/)
* [pandas](https://pandas.pydata.org/getpandas.html)
* [pandas-plink](https://github.com/limix/pandas-plink)
* [tqdm](https://github.com/tqdm/tqdm)

We recommend running S-PCGC via the [Anaconda Python distribution](https://www.anaconda.com/download/). In Anaconda, you can install all the packages with the command "conda install \<package_name\>". Alternatively, the packages can be installed with the command "pip install --user \<package_name\>".

Once all the prerequisite packages are installed, you can install S-PCGC on a git-enabled machine by typing:
```
git clone https://github.com/omerwe/S-PCGC
```

Alternatively, you can download S-PCGC as a zip file from the Github website.
After downloading, we recommend checking that everything is ok by typing ```python test_sumstats.py``` (if you have [pytest](https://github.com/pytest-dev/pytest) installed, you can simply go the project directory and type ```pytest```). It should run for a few minutes. If everything is ok then you shouldn't see an error.


<br><br>
# Usage overview
S-PCGC performs a case-control analysis in four stages:
1. **Generate a sync file for your annotations**. This is a very simple offline step that only needs to be run once. It only gathers some information about the annotations (e.g., the minimum value of each annotation across all SNPs).
2. **Estimate the cross-product of r^2 values across all pairs of functional annotations, via a reference panel such as 1000 genomes**. This step is similar to LD-score computation.
3. **Generate summary statistics**. These are specialized summary statistics explicitly designed for S-PCGC (unlike standard summary statistics analyzed by S-LDSC).
4. **Estimate functional enrichment**. This is what we came here for.

S-PCGC fully supports S-LDSC input and output formats, which enables it to interact with the S-LDSC ecosystem. We therefore recommend that you [familiarize yourself with S-LDSC](https://github.com/bulik/ldsc/wiki) before running S-PCGC.
For example, you can input S-PCGC results into [the S-LDSC scripts for partitioned heritability from continuous annotations](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability-from-Continuous-Annotations). 

<br><br>
# A toy example
The following is a simple end-to-end S-PCGC analysis, which can run in ~2 minutes. We will estimate heritability, genetic correlation and enrichment for four simulated functional annotations in two simulated case-control studies with 250 shared controls and a disease population prevalence of 1%. To estimate the cross-product of r^2 values, we will  use a (simulated) reference panel. All the input files are found in the directory `example`. To begin, please cd into the S-PCGC directory, and type the following commands (using the anaconda version of python if available):
```
mkdir temp_results

#create a sync file (a one-time offline operation)
python pcgc_sync.py --annot-chr example/model. --out temp_results/model

#Compute cross-r^2 between functional annotations
python pcgc_r2.py \
--bfile example/ref_panel \
--annot-chr example/model. \
--sync temp_results/model. \
--out temp_results/prodr2

#Compute summary statistics for study 1
python pcgc_sumstats_creator.py \
--bfile example/s1 \
--pheno example/s1.phe \
--covar example/s1.cov \
--frqfile-chr example/model. \
--annot-chr example/model. \
--sync temp_results/model. \
--prev 0.01 \
--out temp_results/s1

#Compute summary statistics for study 2
python pcgc_sumstats_creator.py \
--bfile example/s2 \
--pheno example/s2.phe \
--covar example/s2.cov \
--frqfile-chr example/model. \
--annot-chr example/model. \
--sync temp_results/model. \
--prev 0.01 \
--out temp_results/s2

#Run S-PCGC to estimate h^2, rg and functional enrichment
python pcgc_main.py \
--annot-chr example/model. \
--sync temp_results/model. \
--frqfile-chr example/model. \
--sumstats temp_results/s1.,temp_results/s2. \
--prodr2 temp_results/prodr2. \
--out temp_results/results

#view heritability estimates for the two studies
cat temp_results/results.s1.output | column -t
cat temp_results/results.s2.output | column -t

#view a table of genetic correlation estimates between the studies
cat temp_results/results.rg

#view the functional enrichment estimates of the two studies
cat temp_results/results.s1.results | column -t
cat temp_results/results.s2.results | column -t
```

#### Some quick comments about this example:
1. Many of the flags are analogous to S-LDSC flags and have similar names.
2. The `.results` files have the exact same format as S-LDSC output files.
3. The `.output` results show both marginal and conditional heritability (please see details below).
4. The flag ```--sumstats``` can accept any number of comma-separated files.
5. In this example we ran `pcgc_r2.py` and `pcgc_sumstats_creator.py` on the entire genome. In real data analysis it may be easier to run these scripts on each chromosome separately, and then use the flag `--prodr2-chr`, `--sumstats-chr` when calling `pcgc_main.py`
6. S-PCGC supports many more options than shown here. For a full list and explanations, please type ```python <file_name> --help```


# An example with real annotations
The following example uses simulated genotypes and real functional annotations from the [Baseline-LD model](https://www.nature.com/articles/ng.3954). In this example, we will use a 'representative set of genotyped or well-imputed SNPs' stored in the file example/good_snps.txt. (generally, 'good SNPs' should be genotyped or well-imputed SNPs (e.g. INFO score >0.9) that passed QC and are not in the MHC region --- see details below).
To run this example, you need a directory called `1000G` that contains plink files of European individuals from the 1000 genomes reference panel (one file per chromosome; see download instructions below). You can create a symbolic link to this directory from your working directory with the command `ln -s <path_to_1000G> 1000G`.
```
#download and uncompress the Baseline-LD (v2) annotations
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v2.0_ldscores.tgz
tar -xzvf 1000G_Phase3_baselineLD_v2.0_ldscores.tgz

#run pcgc_sync.py to collect annotations details
python pcgc_sync.py --annot-chr baselineLD_v2.0/baselineLD. --out baselineLD_v2.0/baselineLD

#run pcgc_r2.py on each chromosome file, using the set of 'good SNPs'
for i in {1..22};
do
    python pcgc_r2.py \
    --annot baselineLD_v2.0/baselineLD.${i}. \
    --sync baselineLD_v2.0/baselineLD. \
    --bfile 1000G/1000G.EUR.QC.${i} \
    --extract example/good_snps.txt \
    --out baselineLD_v2.0/baselineLD.goodSNPs.${i} 
done

#Compute 1000G MAFs (please change ~/plink/plink to the local path of your plink executable)
for i in {1..22};
do
    ~/plink/plink \
    --bfile 1000G/1000G.EUR.QC.${i} \
    --freq \
    --out 1000G/1000G.EUR.QC.${i}
done

#Create summary statistics
mkdir s1_sumstats
for i in {1..22};
do
    python pcgc_sumstats_creator.py \
    --bfile example/s1_chr${i} \
    --extract example/good_snps.txt \
    --pheno example/s1.phe \
    --covar example/s1.cov \
    --annot baselineLD_v2.0/baselineLD.${i}. \
    --sync baselineLD_v2.0/baselineLD. \
    --prev 0.01 \
    --frqfile 1000G/1000G.EUR.QC.${i}. \
    --out s1_sumstats/s1_chr${i}
done

#estimate heritability and functional enrichment
python pcgc_main.py \
--annot-chr baselineLD_v2.0/baselineLD. \
--sync baselineLD_v2.0/baselineLD. \
--frqfile-chr 1000G/1000G.EUR.QC. \
--sumstats-chr s1_sumstats/s1_chr \
--prodr2-chr baselineLD_v2.0/baselineLD.goodSNPs. \
--out s1_sumstats/pcgc

#view heritabiltiy estimates
cat s1_sumstats/pcgc.s1_chr.output | column -t

#view functional enrichment estimates
cat s1_sumstats/pcgc.s1_chr.results | column -t

```


<br><br>
# Additional information

<br><br>
## Obtaining reference panel and annotation files

S-PCGC is fully compatible with the S-LDSC input format. It requires two pieces of information also used by S-LDSC:
1. Functional annotation files. We recommend using the Baseline-LD 2.0 model [(Gazal et al. 2017)](https://www.nature.com/articles/ng.3954). You can download it by typing:
```
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v2.0_ldscores.tgz
```
2. A reference panel (required to compute cross-product of r^2 values). We recommend using 1000 genomes data from a relevant population. For example, you can download 1000G data from [the 1000 Genomes FTP site](https://bit.ly/2OyfNaL) and then convert it to plink format using the [`plink --vcf`](https://www.cog-genomics.org/plink2/input#vcf) command.

<br><br>
## Regressing genotypes on principal components
It is common to include principal components (PCs) as covariates to prevent possible confounding due to population structure. We recommend computing PCs via external software (e.g. [FlashPCA2](https://github.com/gabraham/flashpca)) and including them as additional covariates before computing summary statistics.

A complexity of case-control studies is that the PCs cannot be regressed from the phenotypes. We therefore recommend to regress PCs from the genotypes. You can do this in `pcgc_sumstast_creator.py` via the following two flags:
1. `--covars-regress <covariate names>`. This is a comma-separated list of names of covariates that will be regressed out of the genotypes (i.e., PCs).
2. `--pve <pve file>`. This is the name of a file with a single column and no header that shows the fraction of variance explained by each PC mentioned in `--covars-regress`. Common PCA software like [FlashPCA2](https://github.com/gabraham/flashpca) produce such a file.

<br><br>
## Case-Control GWAS Simulator
To gain more intuition about case-control studies, you may wish to look at the GWAS simulator code included in this package. As a gentle interface, you can look at the last few lines of the file `run_pcgc_simulations.py`. These lines continuously simulate case-control studies with annotations and with no LD, apply PCGC and GCTA estimates, and report the average estimation bias compared to the true generative values. S-LDSC is not included in the comparison because it cannot be invoked in the absence of LD (simulations with LD are possible but are much slower and more complicated; please contact me if you are interested in these).


<br><br>
# Important notes
1. `pcgc_r2.py` must use **exactly** the same SNPs as `pcgc_sumstats_creator.py`. If you use `--extract` in one of them, you must use it in the other as well.

2. The set of SNPs with summary statistics (as in `good_snps.txt` in the example above) should ideally be a set of genotyped or well-imputed SNPs (e.g. INFO score > 0.9) that passed QC. It is highly recommended to also exclude the MHC region (chromosome 6 28M-32M) from these SNPs. One reasonable choice of SNPs is the set of HapMap3 SNPs (which you can download [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2)), as these SNPs are typically genotyped or well-imputed. However, please make sure to also exclude SNPs within the MHC region.

3. Overlapping individuals (shared between the two studies) will not be automatically detected. Please make sure that overlapping individuals are clearly marked in the plink files by having exactly the same family id and individual id.

4. To account for population stratification, it is not enough to simple include principal components as covariates as is typically done. Instead, you need to invoke ```pcgc_sumstats_creator.py``` with the flags ```--covars-regress``` and ```--pve``` (please see the detailed explanation above).

5. The code assumes that the first annotation is always a base annotation that simply assigns 1.0 to all SNPs. Please adhere to this convention!

6. The .output files report two numbers: marginal and conditional heritability. The difference is that conditional heritability conditions on the covariates in the study, meaning that it ignores the variance introduced by the covariates. Most studies report conditional heritability, but we believe that the marginal heritability is more informative. For details, please see the [PCGC-s paper](https://www.sciencedirect.com/science/article/pii/S0002929718301952).

7. We highly recommend that you provide external estimates of SNP frequencies to pcgc_sumstats_creator.py via the flag `--frqfile`, based on a reference panel. This can prevent bias arising due to the fact that cases are over-enriched in case-control studies, which could bias the MAF estimates and the heritability estimates.



<br><br>
# FAQ
Q: Can I create my own annotations?<br>
A: Yes! Please read the section called "Creating an annot file" in the [S-LDSC wiki](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial) for instructions.

Q: Can I run just a basic analysis without functional annotations?<br>
A: Yes. Just omit the --annot and --annot-chr flags, and S-LDSC will perform a simple analysis.

Q: Is it safe to make the summary statistics files publicly available?<br>
A: Absolutely! These summary statistics don't expose any individual-level data. You may see that some of the output files include individual ids, but this is only intended to help identify overlapping individuals. The only information included about these individuals is (a) the noise in the estimation of their kinship with themselves (which should be 1.0 in expectation), and (b) the noise multiplied by their phenotype value. These fields are required to obtain accurate heritability/rg estimates (please see the [PCGC-s paper](https://www.sciencedirect.com/science/article/pii/S0002929718301952) for details)

Q: Should I include imputed SNPs in the summary statistics?<br>
A: Not necessarily. The set of SNPs with summary statistics should be a representative set of common SNPs, correspond to the "regression SNPs" of S-LDSC. We recommend using a set of reliable common SNPs that are either genotyped or well-imputed, such as HapMap3 SNPs (which you can download [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2)). We also recommend excluding SNPs within the MHC (chromosome 6 28M-32M) from all analyses, as done by S-LDSC and other tools.

Q: Can I use my case-control data for `pcgc_r2.py` instead of a reference panel?<br>
A: Yes, but this has two caveats: (1) the analysis will be slower with larger data, and (2) for genetic correlation, it makes more sense to use an external reference panel instead of arbitrarily choosing one of the two studies. If you only care about heritability, there's no need for a reference panel. In this case, we suggest that you downsample your data to ~5000 individuals and only run `pcgc_r2.py` on this subset for greater speed.

Q: Can I use standard (publicly available) summary statistics instead of having to create my own summary statistics?<br>
A: Unfortunately no. To obtain unbiased estimates, S-PCGC must uses specialized summary statistics.

Q: Can my data include related individuals?<br>
A: No.

Q: Does S-PCGC recognize categorical covariates?<br>
A: No (though we might add this in the future). Please transform all your categorical covariates into a series of binary dummy variables.

Q: Can `pcgc_sumstats_creator.py` use bgen files instead of plink files?<br>
A: Currently no.

Q: Can my data include two studies with overlapping individuals?<br>
A: Yes, as long as such individuals are clearly marked in the plink files by having exactly the same family id and individual id. Otherwise, you might get severely biased results.

Q: Can I compute the rg for each annotation separately?<br>
A: Yes, by using the flag `--rg-annot`. This will create a separate .rg file for every pair of studies. However, please note that the results may be nonsensical for continuous annotations, because they can explain a negative amount of heritability.

Q: Can S-PCGC estimate heritability directly from raw genotypes, without using summary statistics?<br>
A: No. In our experience using summary statistics is preferable, because it allows extremely fast performance at a negligible loss of accuracy. However, if you want an exact PCGC implementation, we recommend trying out [LDAK](http://dougspeed.com/pcgc-regression/). Note that the LDAK implementation is limited to less than 100,000 individuals and 20 annotations.

Q: Can S-PCGC fit an intercept like LDSC?<br>
A: No. There's no need to estimate an intercept when the summary statistics are created with `pcgc_sumstats_creator.py`, because the intercept is already known.

<br><br>
-----------------
Contact
---------
For questions and comments, please contact Omer Weissbrod at oweissbrod[at]hsph.harvard.edu



