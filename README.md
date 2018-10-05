# S-PCGC
Heritability, genetic correlation and functional enrichment estimation for case-control studies

S-PCGC is an adaptation of stratified LD score regression [(S-LDSC)](https://www.nature.com/articles/ng.3404) for case-control studies. Similar to S-LDSC, S-PCGC can estimate genetic heritability, genetic correlation and functional enrichment. However, S-PCGC is explicitly designed for case-control studies rather than quantitative phenotypes. This can make a large differnce, especially in the presence of strong non-genetic risk factors. Parts of the S-PCGC code are adapted from S-LDSC with permission.

The main features of S-PCGC are:
1. **S-PCGC is designed for case-control studies**. Such studies include many subtleties not accounted for by methods designed for quantitative phenotypes like S-LDSC.
2. **Seamless integration with S-LDSC format**. S-PCGC accepts the same input files and outputs the same output files as S-LDSC.
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
# Usage Overview
S-PCGC carries performs a case-control analysis in four stages:
1. **Generate a sync file for your annotations**. This is a very simple offline step that only needs to be run once. It only gathers some information about the annotations (e.g., the minimum value of each annotation across all SNPs).
2. **Generate summary statistics**. These are specialized summary statistics explicitly designed for S-PCGC (unlike standard summary statistics analyzed by S-LDSC).
3. **Estimate the cross-product of r^2 values across all pairs of functional annotations, via a reference panel such as 1000 genomes**. This step is similar to LD-score computation. Note that this step requires a summary statistics file *only* to identify which SNPs have summary statistics --- the actual summary statistics are not used.
4. **Estimate functional enrichment**. This is what we came here for.

S-PCGC fully supports S-LDSC input and output formats, which enables it to interact with the S-LDSC ecosystem. We therefore recommend that you [familiarize yourself with S-LDSC](https://github.com/bulik/ldsc/wiki) before running S-PCGC.
For example, you can input S-PCGC results into [the S-LDSC scripts for partitioned heritability from continuous annotations](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability-from-Continuous-Annotations). 

<br><br>
## TL;DR (A Simple Example)
The following is a simple end-to-end S-PCGC analysis, which can run in ~2 minutes. We will estimate heritability, genetic correlation and functional enrichment for two simulated case-control studies with 250 shared controls and a disease population prevalence of 1%, using only SNPs in chromosome 1. To estimate the cross-product of r^2 values, we will  use a (simulated) reference panel. All the input files are found in the directory `example`. To begin, please cd into the S-PCGC directory, and type the following commands (using the anaconda version of python if available):
```
mkdir temp_results

#create a sync file (a one-time offline operation)
python pcgc_sync.py --annot-chr example/model. --out temp_results/model

#Compute summary statistics for study 1
python pcgc_sumstats_creator.py \
--bfile example/s1 \
--pheno example/s1.phe \
--covar example/s1.cov \
--frqfile example/model.1. \
--annot example/model.1. \
--sync temp_results/model. \
--prev 0.01 \
--out temp_results/s1

#Compute summary statistics for study 2
python pcgc_sumstats_creator.py \
--bfile example/s2 \
--pheno example/s2.phe \
--covar example/s2.cov \
--frqfile example/model.1. \
--annot example/model.1. \
--sync temp_results/model. \
--prev 0.01 \
--out temp_results/s2

#Compute cross-r^2 between functional annotations, using only SNPs with summary statistics
python pcgc_r2.py \
--bfile example/ref_panel \
--annot example/model.1. \
--sync temp_results/model. \
--sumstats temp_results/s1.sumstats.gz \
--out temp_results/prodr2.1

#Run S-PCGC to estimate h^2, rg and functional enrichment
python pcgc_main.py \
--annot example/model.1. \
--sync temp_results/model. \
--frqfile example/model.1. \
--sumstats temp_results/s1.,temp_results/s2. \
--prodr2 temp_results/prodr2.1. \
--out temp_results/results

#view heritability estimates for the two studies
cat temp_results/results.s1.output
cat temp_results/results.s2.output

#view a table of genetic correlation estimates between the studies
cat temp_results/results.rg

#view the functional enrichment estimates of the two studies
cat temp_results/results.s1.results
cat temp_results/results.s2.results
```

#### Some quick comments about this example:
1. Many of the flags are analogous to S-LDSC flags and have similar names.
2. The `.results` files have the exact same format as S-LDSC output files.
3. The `.output` results show both marginal and conditional heritability (please see details below).
4. The flag ```--sumstats``` can accept any number of comma-separated files.
4. If we wanted to run a whole-genome analysis instead of just chromosome 1, we just need to replace the flags `--annot`, `--frqfile` with `--annot-chr`, `--frqfile-chr`, and remove the `.1` suffix from all the input and output files (check it out!)
5. We can also replace the flags `--annot`, `--frqfile` with `--annot-chr`, `--frqfile-chr` only in the last command (`pcgc_main.py`). This will cause S-PCGC to use *only* chromosome 1 SNPs to estimate annotation effects, but to report heritability and enrichment using all the common SNPs found in the reference panel.
6. S-PCGC supports many more options than shown here. For a full list and explanations, please type ```python <file_name> --help```



<br><br>
# Obtaining Reference Panel and Annotation files

S-PCGC is fully compatible with the S-LDSC input format. It requires two pieces of information also used by S-LDSC:
1. Functional annotation files. We recommend using the Baseline-LD 2.0 model [(Gazal et al. 2017)](https://www.nature.com/articles/ng.3954). You can download it by typing:
```
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v2.0_ldscores.tgz
```
2. A reference panel (required to compute cross-product of r^2 values). We recommend using 1000 genomes data from a relevant population. For example, you can download 1000G data from [the 1000 Genomes FTP site](https://bit.ly/2OyfNaL) and then convert it to plink format using the [`plink --vcf`](https://www.cog-genomics.org/plink2/input#vcf) command.

<br><br>
## Regression of Genotypes on Principal Components
It is common to include principal components (PCs) as covariates in the analysis to prevent possible confounding due to population structure. We recommend computing PCs via external software (e.g. [FlashPCA2](https://github.com/gabraham/flashpca)) and including them as additional covariates before computing summary statistics.

A complexity of case-control studies is that the PCs cannot be regressed from the phenotypes. We therefore recommend to regress PCs from the genotypes. This can be done in `pcgc_sumstast_creator.py` via the following two flags:
1. `--covars-regress <covariate names>`. This should be a comma-separated list of names of covariates that be regressed on the genotypes (i.e., PCs).
2. `--pve <pve file>`. This should be the name of a file with a single column and no header that shows the fraction of variance explained by each PC mentioned in `--covars-regress`. Common PCA software like [FlashPCA2](https://github.com/gabraham/flashpca) produce such a file.


<br><br>
# Important notes
1. Overlapping individuals (shared between the two studies) will not be automatically detected. Please make sure that overlapping individuals are clearly marked in the plink files by having exactly the same family id and individual id.

2. To account for population stratification, it is not enough to simple include principal components as covariates as is typically done. Instead, you need to invoke ```pcgc_sumstats_creator.py``` with the flags ```--covars-regress``` and ```--pve``` (please see the detailed explanation above).

3. The code assumes that the first annotation is always a base annotation that simply assigns 1.0 to all SNPs. Please adhere to this convention!

4. The .output files report two numbers: marginal and conditional heritability. The difference is that conditional heritability conditions on the covariates in the study, meaning that it ignores the variance introduced by the covariates. Most studies report conditional heritability, but we believe that the marginal heritability is more informative. For details, please see the [PCGC-s paper](https://www.sciencedirect.com/science/article/pii/S0002929718301952).

5. We highly recommend that you provide external estimates of SNP frequencies to pcgc_sumstats_creator.py via the flag `--frqfile`, based on a reference panel. This can prevent bias arising due to the fact that cases are over-enriched in case-control studies, which could bias the MAF estimates and the heritability estimates.

6. S-PCGC estimates h^2 across all common SNPs in the reference panel, including ones without summary statistics (e.g. unimputed SNPs), exactly like S-LDSC. This stands in contrast to software such as [gcta](https://cnsgenomics.com/software/gcta/), which only estimates the heritability that is tagged by genotyped and imputed SNPs.



<br><br>
# FAQ
Q: Can I create my own annotations?<br>
A: Yes! Please read the section called "Creating an annot file" in the [S-LDSC wiki](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial) for instructions.

Q: Can I run just a basic analysis without functional annotations?<br>
A: Sure! Just omit the --annot and --annot-chr flags, and S-LDSC will perform a simple analysis.

Q: Is it safe to make the summary statistics files publicly available?<br>
A: Absolutely! These summary statistics don't expose any individual-level data. You may see that some of the output files include individual ids, but this is only intended to help identify overlapping individuals. The only information included about these individuals is (a) the noise in the estimation of their kinship with themselves (which should be 1.0 in expectation), and (b) the noise multiplied by their phenotype value. These fields are required to obtain accurate heritability/rg estimates (please see the [PCGC-s paper](https://www.sciencedirect.com/science/article/pii/S0002929718301952) for details)

Q: Should I include imputed SNPs in the summary statistics?<br>
A: Not really. S-PCGC first estimates annotation effects via a subset of SNPs, and then estimates heritability and enrichment for all common SNPs that appear in the reference panel. The subset of SNPs corresponds to the set of "regression SNPs" of S-LDSC. We recommend using a set of high-quality common SNPs, such as HapMap3 SNPs. We also recommend excluding SNPs within the MHC (chromosome 6 28M-32M) from all analyses, as done by S-LDSC and other software tools.

Q: Can I use standard (publicly available) summary statistics instead of having to create my own summary statistics?<br>
A: Unfortunately no. To obtain unbiased estimates, S-PCGC must uses specialized summary statistics.

Q: Can my data include related individuals?<br>
A: No.

Q: Does S-PCGC recognize categorical covariates?<br>
A: No (though we might add this in the future). Please transform all your categorical covariates into a series of binary dummy variables.

Q: Can pcgc_sumstats_createor use bgen files instead of plink files?<br>
A: Currently no.

Q: Can my data include two studies with overlapping individuals?<br>
A: Yes, as long as such individuals are clearly marked in the plink files by having exactly the same family id and individual id. Otherwise, you might get severely biased results.

Q: Can I compute the rg for each annotation separately?<br>
A: Yes, by using the flag `--rg-annot`. This will create a separate .rg file for every pair of studies. However, please note that the results may be nonsensical for continuous annotations, because they can explain a negative amount of heritability.

Q: Can S-PCGC estimate heritability without using summary statistics?<br>
A: No. In our experience using summary statistics is preferable, because it allows extremely fast performance at a negligible loss of accuracy. However, if you want an exact PCGC implementation, we recommend trying out [LDAK](http://dougspeed.com/pcgc-regression/). Note that the LDAK implementation is limited to less than 100,000 individuals and 20 annotations.



<br><br>
-----------------
Contact
---------
For questions and comments, please contact Omer Weissbrod at oweissbrod[at]hsph.harvard.edu



