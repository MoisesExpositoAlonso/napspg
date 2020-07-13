# Supplementary Notes on Non-Additive Polygenic Models (NAPs)
by Moi Exposito-Alonso & Rasmus Nielsen. Last updated: Jul 13 2020

#### Note 1: The interpretation of fitness, its natural scale, and typical datasets
Darwinian fitness includes both viability and fecundity components. Fitness thus has a natural scale bounded at zero (individuals that did not survive or did not produce offspring) to any positive integer number.  Such data behaves as a Poisson or negative binomial in species where adults produce very few offspring, and approximates a Normal in species that have large numbers of offspring. The species analyzed in our original publication, yeast, fruit fly, and thale cress, produce high numbers and thus we chose to model fitness of a population as a combination of Normal distributions with genotype-specific mean and variance variances, and added a zero-inflation term to account for the many zeroes, common in fitness datasets. This approach would likely be more appropriate to model other traits than fitness that are also bounded to zero and with heavy tails (number of morphological structures, etc.). However, we must note that when a trait does not have a natural measurement scale (certain concentrations of metabolites, or indirect absorbance metrics), global epistatic models may fit best the data but may only indicate a non-linearity in the measuring scale.

#### Note 2: Likelihood optimization and file IO
We re-implemented the Spatial Projected Gradient method [1] from the original source C++ code by the authors (available at www.ime.usp.br/~egbirgin/tango), and added a Spectral steplength subroutine that improves convergence [3]. In order to read large genome matrices fast, we refactored the C++ code from plink2R R package (https://github.com/gabraham/plink2R) to read PLINK’s .bed binary files into a C++ Arma SNP matrix of 0, 1, and 2, for reference homozygote, heterozygote, and alternative homozygote. Missing values were treated as reference alleles, as it is expected that a selection of high-quality and strongly significant SNPs are used as input.

#### Note 3: The use of individual cross-validation for model selection
We decided to focus on evaluating models based on cross-validation accuracy instead of likelihood because cross-validation accuracy is commonly used in evaluating genomic prediction and polygenic risk scores. We must also note that because fitness is bounded at zero, and the additive model—being poorly defined—makes negative fitness predictions, the likelihood under a multiplicative model typically outperformed additive models even when the multiplicative model underperformed the additive model in cross-validation accuracy. In addition, optimization issues that differ between additive and multiplicative models could bias our conclusions if we had focused on likelihood maximization. Future development testing a larger variety of global epistasis models and further likelihood optimization will likely be fruitful to understand the poor likelihood of additive models.

#### Note 4: Global epistasis parameter inference
Our optimization algorithm has the ability to optimize the exact value of the global epistatic hyperparameter. We, however, opted for a grid search approach of 9 values (from 0.8 to 1.2 in 0.05 increments) to speed up optimization. 

#### Note 5: Selection coefficients inference
While our optimization returns a list of selection coefficients of each of the input SNPs, we focus on testing the predictability at the individual level because it is at such level that global epistasis patterns become apparent. During our method development, we nevertheless used the inferences of selection coefficients in simulated datasets to confirm that our optimization correctly converged to the simulated selection coefficients. The inference of selection coefficients, and how it is affected by global epistasis and other varying SNPs, is something we aim to continue studying and improving in the future.

#### Note 6: On multi-locus association with BSLMM
The BSLMM GWA method implemented in GEMMA [3] is a necessary pre-step, as efficient implementations of linear models allow to pre-screen millions of SNPs. Furthermore, it ameliorates the potential problem of pre-selecting SNPs in high LD. Because BSLMM fits all SNPs at once, the effects of two alleles in linkage disequilibrium should ideally not be double-estimated, but the one with the higher marginal effect will be favored as such marginal effect is used as prior. If they are 100% in LD, the effect of one SNP is set to 0 at random. Alternatively, or in addition, one could pre-select the leading SNP from GWA for each haploblock in the genome. While well-defined haploblocks exist in the human genome [4], but out-of-the box definitions are not available in other species, and we preferred to avoid an additional pre-step to run NAP.

#### References
[1] Birgin, E., J. Martínez, and M. Raydan. 2000. “Nonmonotone Spectral Projected Gradient Methods on Convex Sets.” SIAM Journal on Optimization: A Publication of the Society for Industrial and Applied Mathematics 10 (4): 1196–1211. https://doi.org/10.1137/S1052623497330963. \
[2] R Varadhan and C Roland (2008), Simple and globally-convergent methods for accelerating the convergence of any EM algorithm, Scandinavian J Statistics, doi: 10.1111/j.1467-9469.2007.00585.x. \
[3] Zhou, Xiang, Peter Carbonetto, and Matthew Stephens. 2013. “Polygenic Modeling with Bayesian Sparse Linear Mixed Models.” PLoS Genetics 9 (2): e1003264. https://doi.org/10.1371/journal.pgen.1003264. \
[4] Berisa, Tomaz, and Joseph K. Pickrell. 2016. “Approximately Independent Linkage Disequilibrium Blocks in Human Populations.” Bioinformatics 32 (2): 283–85. https://doi.org/10.1093/bioinformatics/btv546. \



