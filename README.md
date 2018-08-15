# BPSC
**B**eta-**P**oisson model for **S**ingle-**C**ell RNA-seq data analyses

# How to install "BPSC"
### Latest release
[Version 0.99.2](https://github.com/nghiavtr/BPSC/releases/download/v0.99.2/BPSC_0.99.2.tar.gz)

What's new in version 0.99.2
- Use t-test for weird isoform data, for example all samples of the control group are unexpressed, to avoid NA results in BPglm. Great thanks to Tian Mou for the contribution!
- Use RMarkdown instead of Sweave for vignette document.

What's new in version 0.99.1
- Fix and improve vignette documents
- Use quasi-Possion as the family function if the glm fitting (in function BPglm) with BPfam is not converged

#### Older versions can be downloaded from here https://github.com/nghiavtr/BPSC/releases
- [Version 0.99.1](https://github.com/nghiavtr/BPSC/releases/download/v0.99.1/BPSC_0.99.1.tar.gz)
- [Version 0.99.0](https://github.com/nghiavtr/BPSC/releases/download/v0.99.0/BPSC_0.99.0.tar.gz)

##### Install from command line:
```R
R CMD INSTALL BPSC_x.y.z.tar.gz 
```
where BPSC_x.y.z.tar.gz is one version of BPSC
##### BPSC package requires the some dependent packages:
```R
statmod, doParallel
```
### Development version from Github using R:
```R
install.packages("devtools")
library("devtools")
install_github("BPSC","nghiavtr")
library("BPSC")
```
# View vignette for user guide
```R
vignette("BPSC")
```
# Summary

Single-cell RNA-sequencing technology allows detection of gene expression at the single-cell level. One typical feature of the data is a bimodality in the cellular distribution even for highly expressed genes, primarily caused by a proportion of non-expressing cells. The standard and the over-dispersed gamma-Poisson models that are commonly used in bulk-cell RNA-sequencing are not able to capture this property.

We introduce a beta-Poisson mixture model that can capture the bimodality of the single-cell gene expression distribution. We further integrate the model into the generalized linear model framework in order to perform differential expression analyses. The whole analytical procedure is called BPSC. The results from several real single-cell RNA-seq datasets indicate that ~90% of the transcripts are well characterized by the beta-Poisson model; the model-fit from BPSC is better than the fit of the standard gamma-Poisson model in >80% of the transcripts. Moreover, in differential expression analyses of simulated and real datasets, BPSC performs well against edgeR, a conventional method widely used in bulk-cell RNA-sequencing data, and against scde and MAST, two recent methods specifically designed for single-cell RNA-seq data.

Related publication:
Vu,T.N. et al. (2016) Beta-Poisson model for single-cell RNA-seq data analyses. Bioinformatics, btw202.
http://bioinformatics.oxfordjournals.org/content/early/2016/04/18/bioinformatics.btw202
