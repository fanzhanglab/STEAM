
# STEAM

<img width="30%" align = "right" src="https://github.com/fanzhanglab/STEAM/blob/main/media/STEAM_logo.png?raw=true">

<!-- badges: start -->

[![R-CMD-check](https://github.com/fanzhanglab/STEAM/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/fanzhanglab/STEAM/actions/workflows/check-standard.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Visitors](https://api.visitorbadge.io/api/visitors?path=https%3A%2F%2Fgithub.com%2Ffanzhanglab%2FSTEAM-v4&label=%23Visits&labelColor=%23000000&countColor=%2300c00B&style=plastic)

<!-- badges: end -->

<p align="justify">

One common challenge in evaluating the robustness of identified cell
type clusters in spatial omics is the lack of ground truth cell type
labels in real-world data from disease conditions. To address this, we
introduce STEAM, a Spatial Transcriptomics Evaluation Algorithm and
Metric for clustering performance, developed to evaluate the consistency
of clustering algorithms and the reliability of cell annotations in
spatial omics data. Our hypothesis is that if clusters are robust and
consistent across tissue regions, selecting a subset of cells or spots
within a cluster should enable accurate prediction of cell type
annotations for the remaining cells within that cluster, due to spatial
proximity and gene expression covarying patterns.

STEAM incorporates various machine learning models, including Random
Forest, XGBoost, and SVM, to assess the prediction accuracy and
consistency of clusters, along with statistical metrics like Kappa
score, F1 score, ARI (Adjusted Rand Index), etc. We demonstrated the
capability of STEAM on multi-cell and single-cell resolution spatial
transcriptomics and proteomics. Notably, STEAM supports multi-sample
training, enabling the evaluation of cross-replicate clustering
consistency. Furthermore, we used STEAM to evaluate the performance of
spatial-aware and spatial-ignorant clustering methods, offering
researchers a valuable tool for more informed result interpretation.
</p>

<img width="100%" align = "center" src="https://github.com/fanzhanglab/STEAM/blob/main/media/Figure1.png?raw=true">

<!-- <img width="100%" align = "center" src="https://github.com/fanzhanglab/STEAM/blob/main/man/figures/Figure1.png"> -->

</br>

## Installation

You can install the STEAM Package from
[GitHub](https://github.com/fanzhanglab/STEAM/) using the devtools as
follows:

``` r
# install.packages("devtools")
devtools::install_github("fanzhanglab/STEAM")
```

(OR)

``` r
remotes::install_github("fanzhanglab/STEAM")
```

<br/>

### Dependencies / Other required packages

``` r
- R (\>= 4.2)
- ggplot2 (\>= 3.4.2)
- caret
- randomForest
- e1071
- scales
- gridExtra
- grid
- reshape2
- viridis
```

<br/>

## Tutorials

**Step-by-step notebook** of applying STEAM on 10X Visium Human Brain
Data (DLPFC):

- <a href="https://htmlpreview.github.io/?https://github.com/fanzhanglab/STEAM/blob/main/tutorials/STEAM_DLPFC_tutorial.html">
  Tutorial of applying STEAM on DLPFC data </a>

#### Below are several major steps of running STEAM:

``` r
# Create a new STEAM object for the loaded spatial transcriptomic data
STEAM.Obj <- LoadSTEAM(count_exp = matrix, spatial = coordinates, labels = labels, Seurat.obj = NULL)
```

``` r
STEAM.Obj <- RunSTEAM(STEAM.obj, train.ratio = 0.8, n.size = 5, seed = 123, cv.folds = 10, cv.repeats = 3, trainval.ratio = 0.8, model = "rf", n.tree = 500, kernel = 'linear', train.folder.name = 'train.out', allowParallel = FALSE)
```

<br/>

## Citations

Reynoso, S., Schiebout, C., Krishna, R., Zhang, F. STEAM: Spatial
Transcriptomics Evaluation Algorithm and Metric for clustering
performance, [*bioRxiv*](https://www.biorxiv.org/content/10.1101/2025.02.17.636505v1), 2025, https://doi.org/10.1101/2025.02.17.636505

<br/>

## Help, Suggestion and Contribution

Using github [**issues**](https://github.com/fanzhanglab/STEAM/issues)
section, if you have any question, comments, suggestions, or to report
coding related issues of STEAM is highly encouraged than sending emails.

- Please **check the GitHub
  [issues](https://github.com/fanzhanglab/STEAM/issues)** for similar
  issues that has been reported and resolved. This helps the team to
  focus on adding new features and working on cool projects instead of
  resolving the same issues!
- **Examples** are required when filing a GitHub issue. In certain
  cases, please share your STEAM object and related codes to understand
  the issues.

<br/>

## Contact

Please contact [fanzhanglab@gmail.com](fanzhanglab@gmail.com) for
further questions or protential collaborative opportunities!
