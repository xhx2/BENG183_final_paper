# scImpute
### By Hanxiang Xu, Kazune Fei, Yanpei Wang


---

## Background and Motivation

Single-cell RNA (scRNA-seq) technologies are now emerging as a powerful tool to capture transcriptoe-wide cell-to-cell variability. ScRNA-seq enables the qunatification of intro-population  heterogeneity at a much higher resolution, potentially revealing dynamics in heterogeneous cell populations and complex tissues. 

One important characteristic of scRNA-seq data is the 'dropout' phenomenon, where a gene is observed at a moderate expression level in one cell but undetected in another cell. The frequency of dropouts is protocol dependent, for instance, droplet-based microfluidics protolcols tend to have a much higher dropout rate compared to the Fluidigm C1 platform. Since this characteristic is dependent on protocol and the accuracy of downstream scRNA-seq analyses-such as differential gene expression analysis, cell-type identification, and recontruction of differential trajectories-relies heavily on accurate gene expression measurements, it is crucial for the statistical or computational methods developed for scRNA-seq to take the dropout issue into condisderation and correct these false zero expression values

Existing imputation methods, such as **MAGIC** and **SAVER**, attempt to adddress this problem but have a significant limitation: they alter all gene expression levels, including those unaffected by dropouts. This could potentially introduce new biases into the data and eliminate meaningful biological variation. The scImpute method is proposed to overcome these limitaitons by simultaneously determining which values are affected by dropout events and performing imputation only on those highly probable dropout entries. This approach ensures that techinical variation resulting from scRNA-seq is reduced, while the biological variation in cells is better preserved, thereby avoiding excess bias. 

## Introduction of scImpute

scImpute is a two-stage statistical method built on the philosophy of targeted imputation and cell similarity. 

<img src="figure_1_scImpute.png" width="1000">

**Figure 1. Overview of the scImpute method** scImpute firstly learns each gene's dropout probability in each cell by fitting a mixture model. Next, scImpute imputes the highly probable dropout values in cell $j$ (gene set $A_j$) by borrowing information of the same gene in other similar cells, which are selected based on gene set $B_j$ (not severely affected by dropout events)

### Identification of likely dropouts

scImpute models the expression of each gene within a cell subpopulation using a mixture model composed of two components. The first component is a Gamma distribution used to account for the dropouts, while the second component is a Normal distribution to represent the actual gene expression levels. For each gene $i$, its expression in cell subpopulation $k$ is modeled as a random variable $X_i^{(k)} with density function: 

$f_{X_i}^{(k)}(x)=\lambda_i^{(k)}\mathrm{Gamma}(x;\alpha_i^{(k)})+(1-\lambda_i^{(k)})\mathrm{Normal}(x;\mu_i^{(k)})$, where $\lambda_i^{(k)}$ is gene $i$'s dropout rate in cell subpopulation $k$, $\alpha_i^{(k)}$, $\beta_i^{(k)}$ are the shape and rate parameters of Gamma distribution, and $\mu_i^{(k)}$, $\alpha_i^{(k)}$ are the mean and standard deviation of Normal distribution. The intuition behind this mixture model is that if a gene has high expression and low vairation in most cells within a cell subpopulation, a zero count is more likely to be a dropout value; on the other hand, if a gene has constantly low or medium expression with high variation, then a zero count may reflect real biological variability. The parameters in mixture model are estimated by the Expectation-Maximization (EM) algorithm, and we denote their estimate as $\hat{\lambda_i^{(k)}}$, $\hat{\alpha_i^{(k)}}$, $\hat{\beta_i^{(k)}}$, $\hat{\mu_i^{(k)}}$ and $\hat{\sigma_i^{(k)}}$

### Learning Cell Similarity

For each cell $j$, select a gene set $A_j$ in need of imputation based on genes' dropout probabilities in cell $j$: $A_j= \lbrace i:d_{ij} \ge t \rbrace$, where $t$ is a threshold on dropout probabilities. There is also a gene set $B_j=\lbrace i:d_{ij}<t \rbrace$
To impute the values in set $A_j$, scImpute need to borrow information from similar cells. Cell similarity is learned only using the high-confidence, accurately measured genes in set $B_j#.  



### 




