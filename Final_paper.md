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

### Identification of Likely Dropouts

scImpute models the expression of each gene within a cell subpopulation using a mixture model composed of two components. The first component is a Gamma distribution used to account for the dropouts, while the second component is a Normal distribution to represent the actual gene expression levels. For each gene $i$, its expression in cell subpopulation $k$ is modeled as a random variable $X_i^{(k)}$ with density function: 

$f_{X_i}^{(k)}(x)=\lambda_i^{(k)}\mathrm{Gamma}(x;\alpha_i^{(k)})+(1-\lambda_i^{(k)})\mathrm{Normal}(x;\mu_i^{(k)})$, where $\lambda_i^{(k)}$ is gene $i$'s dropout rate in cell subpopulation $k$, $\alpha_i^{(k)}$, $\beta_i^{(k)}$ are the shape and rate parameters of Gamma distribution, and $\mu_i^{(k)}$, $\alpha_i^{(k)}$ are the mean and standard deviation of Normal distribution. The intuition behind this mixture model is that if a gene has high expression and low vairation in most cells within a cell subpopulation, a zero count is more likely to be a dropout value; on the other hand, if a gene has constantly low or medium expression with high variation, then a zero count may reflect real biological variability. The parameters in mixture model are estimated by the Expectation-Maximization (EM) algorithm, and we denote their estimate as $\hat{\lambda_i^{(k)}}$, $\hat{\alpha_i^{(k)}}$, $\hat{\beta_i^{(k)}}$, $\hat{\mu_i^{(k)}}$ and $\hat{\sigma_i^{(k)}}$

### Learning Cell Similarity

For each cell $j$, select a gene set $A_j$ in need of imputation based on genes' dropout probabilities in cell $j$: $A_j= \lbrace i:d_{ij} \ge t \rbrace$, where $t$ is a threshold on dropout probabilities. There is also a gene set $B_j=\lbrace i:d_{ij}<t \rbrace$ which have accurate gene expression with high confidence and do not need imputation. We learn cells' similarities through the gene set $B_j$, and use the non-negative least squares (NNLS) regression:

$$
\hat{ \beta^{(j)} }=\text{argmin}_{\beta^{(j)}} \lVert X_{B_j,i}-X_{B_j,N_j}\beta^{(j)} \rVert^2_2
$$

$$
\text{subject \ to \ } \beta^{(j)} \ge 0
$$

$N_j$ represents the indices of cells that are candidate neighbors of cell $j$. The response $X_{B_j,j}$ is a vector representing the $B_j$ rows in the $j$-th column of $X$, the design matrix $X_{B_j,N_j}$ is a sub-matrix of $X$ with dimensions $\| B_j \| \times \|N_j\|$, and the coefficients $\beta^{(j)}$ is a vector of length $\| N_j \|$. 

$\hat{ \beta^{(j)} }$ represents the weighted contribution of each similar neighbor cell to the expression profile of cell $j$. NNLS has the property of leadning to sparse estimated coefficient vector $\hat{ \beta^{(j)} }$ whose components may have exact zeros, so NNLS can be used to select similar cells of cell $j$ from its neighbors $N_j$. 

### Imputation of Dropout Values

The estimated coefficients $\hat{\beta^{(j)}}$ from the set $B_j$ are used to impute the expression of genes in the set $A_j$ in cell $j$:

$$
\hat{X}_{ij} =
\begin{cases}
X_{ij}, & i \in B_j \quad \text{(Value retained as it is accurate/true)} \\
X_{i, N_j}\,\hat{\beta}^{(j)}, & i \in A_j \quad \text{(Value imputed using neighbor cells' weighted contributions)}
\end{cases}
$$

## Results

### scImpute recovers gene expression affected by dropouts

#### scImpute recovers the true expression of the ERCC spike-in transcripts

The ERCC spike-ins are synthetic RNA molecules with known concentrations which serve as gold standards of true expression levels, so that the read counts can be compared with the ture expression for accuracy evaluation. The dataset contains 3005 cells from the mouse somatosensory cortex region. Figure 2 shows that after imputation, the median correlation (of the 3005 cells) between 57 transcripts' read counts and their true concentrations increases from 0.92 to 0.95, and the minimum correlation increases from 0.81 to 0.89. 

<img src="figure2_ERCC_correlation.png" width="1000">

**Figure 2. Correlation between ERCC spike-ins' counts and their true concentration**. The	two	distributions	show	the	correlations between	the	ERCC	spike-ins'	log10(count+1)	and	log	10(concentration)	in	the	3005	mouse	cortex	cells	(one	correlation	per	cell	for	raw	counts	or	counts	corrected	by	scImpute).

Figure 3 further shows that the read counts and true concentrations also present a stronger linear relationship in every single cell.

<img src="figure3_ERCC_four_cells.png" width="1000">

**Figure 3. scImpute improves the dropouts in the ERCC RNA scripts** The y-axis and x-axis give the ERCC spike-ins’ log10(count+1) and log10
(concentration) in four randomly selected mouse cortex cells. The imputed data present stronger linear relationships between the true concentrations and the observed counts

### scImpute correctly imputes the dropout values of 892 annotated cell cycle genes in 182 embryonic stem cells (ESCs) that had been staged for cell cycle phases (G1, G2M and S)

These genes are know to modulate the cell cycle and are expected to have non-zero expression during different stages of the cell cycle. Before imputation, 22.5% raw counts of the cell cycle genes are zeros, which are highly likely due to dropouts. The data are normalized by sequencing depths instead of ERCC spike-ins. After imputations, most of the dropout values are corrected, and true dynamics of these genes in the cell cycle are recovered （Figs 4 and 5）.

<img src="figure4_ESC_heatmap.png" width="1000"> 

**Figure 4 Heatmaps	showing	the	log10(count+1)	of	the	892	cell	cycle	genes	before	and	after	imputation.** . Rows	correspond	to	genes,	and	columns	correspond	to	cells.

<img src="figure5_ESC_violin1.png" width="1000"> 

**Figure 5 Expression	of	cell	cycle	genes	before	and	after	imputation.** Violin	plots	showing	the	log10(count+1)	of	the	892	cell	cycle	genes	in	the	three	phases	(G1,	G2M,	and	S).	This	comparison	result	shows	that scImpute	has	successfully	imputed the	dropout	expression	values	of	cell	cycle	genes.

Figure 6 shows that imputed counts also represent the true biological variation in these cell cycles genes. 

<img src="figure6_ESC_violin2.png" width="1000"> 

**Figure 6 Violin plots showing the log10(count+1) of nine cell cycle genes** The expression levels of these genes belong to three phases (G1, G2M, and S). scImpute has corrected the dropout values of cell cycle genes. 



## Advantages and Limitation

An attractive advantage of scImpute is that it can be incorporated into most existing pipelines or downstream analysis of scRNA-seq data, such as normalization, differential expression analysis, clustering and classification. 













