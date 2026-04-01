# Medium Experiments Repository

This repository collects small, article-linked subprojects. Each folder is intended to stand on its own: code, notebook, or script first; focused README second.

The top-level README is organized as a simple index so new subprojects can be added without restructuring the whole file.

## Navigation

- [Repository Structure](#repository-structure)
- [Subprojects](#subprojects)
- [ML Metamorphosis](#ml-metamorphosis)
- [Mean-Target Encoding for Tree Splits](#mean-target-encoding-for-tree-splits)
- [Nonlinear Correlation](#nonlinear-correlation)
- [PRIM Segments](#prim-segments)

## Repository Structure

Each subproject lives in its own directory and should usually include:

- a local `README.md` explaining the article connection, main idea, contents, and run instructions
- the code, notebook, or script used for the experiment
- any subproject-specific dependencies or helper modules

Current subprojects:

- [`metamorphosis/`](./metamorphosis)
- [`mt_encoding/`](./mt_encoding)
- [`nonlinear_correlation/`](./nonlinear_correlation)
- [`prim_segments/`](./prim_segments)

## Subprojects

### ML Metamorphosis

- Folder: [`metamorphosis/`](./metamorphosis)
- README: [metamorphosis/README.md](./metamorphosis/README.md)
- Article:
  ["ML Metamorphosis: Chaining ML Models for Optimized Results"](https://medium.com/data-science/ml-metamorphosis-chaining-ml-models-for-optimized-results-d89d952627a9)
- Summary:
  trains a CNN on a small labeled MNIST subset, uses confident pseudo-labels to enlarge the dataset, and then trains a decision tree on the enriched data.

### Mean-Target Encoding for Tree Splits

- Folder: [`mt_encoding/`](./mt_encoding)
- README: [mt_encoding/README.md](./mt_encoding/README.md)
- Article:
  ["Decision Trees Natively Handle Categorical Data"](https://towardsdatascience.com/decision-trees-natively-handle-categorical-data/)
- Summary:
  compares exhaustive category partition search with ordered mean-target splits for decision tree splitting on categorical features.

### Nonlinear Correlation

- Folder: [`nonlinear_correlation/`](./nonlinear_correlation)
- README: [nonlinear_correlation/README.md](./nonlinear_correlation/README.md)
- Article:
  ["An Undeservedly Forgotten Correlation Coefficient"](https://medium.com/data-science/an-undeservedly-forgotten-correlation-coefficient-86245ccb774c)
- Summary:
  evaluates Pearson correlation, a mutual-information-based coefficient `R`, and Chatterjee's `xi` on synthetic datasets covering nonlinearity, symmetry, consistency, and runtime.

### PRIM Segments

- Folder: [`prim_segments/`](./prim_segments)
- README: [prim_segments/README.md](./prim_segments/README.md)
- Article:
  ["Find Unusual Segments in Your Data with Subgroup Discovery"](https://medium.com/data-science/find-unusual-segments-in-your-data-with-subgroup-discovery-2661a586e60c)
- Summary:
  provides a compact implementation of PRIM for subgroup discovery plus a notebook example on bank churn segmentation.
