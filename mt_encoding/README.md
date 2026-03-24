# Mean-Target Encoding for Tree Splits

This subproject accompanies the article ["Decision Trees Natively Handle Categorical Data"](https://towardsdatascience.com/decision-trees-natively-handle-categorical-data/).

The notebook explores a practical question behind that claim: when a decision tree splits on a categorical feature, do we need to test every possible partition of categories, or can we recover the same best split by first ordering categories with mean target encoding and then checking only ordered cut points?

## What is in this folder

- `mt_encoding.ipynb`: a self-contained experiment notebook

## Main idea

For a categorical feature with `k` distinct values, an exhaustive search over binary splits becomes expensive very quickly because the number of possible partitions grows exponentially.

The notebook compares two strategies:

- `all_category_splits`: enumerate all unique binary partitions of the categories
- `ordered_splits`: compute a mean target encoding for each category, sort the categories by that value, and consider only threshold-style splits in that order

It then checks whether both strategies reach the same best split score while measuring the runtime difference.

## What the notebook implements

The notebook contains:

- a synthetic data generator for categorical features with either binary or continuous targets
- split generators for exhaustive search and ordered mean-target splits
- scoring functions for the two settings:
  `gini_impurity` for binary classification and `mse_criterion` for regression
- an experiment loop that repeats each setup 100 times for multiple cardinalities
- plots of average execution time on a log scale

## Experiment design

Each synthetic dataset contains:

- `k` categories
- `100` rows per category
- a target column named `Value`

Two target types are tested:

1. Binary target:
   uses Gini impurity and evaluates `k` from `4` to `12`
2. Continuous target:
   uses mean squared error and evaluates `k` from `4` to `9`

For each `k`, the notebook runs 100 random trials and records:

- whether the best exhaustive score and the best ordered-split score are equal
- average runtime for exhaustive search
- average runtime for ordered splits

## What the results show

The purpose of the experiment is not just speed benchmarking. It checks whether the ordered mean-target approach recovers the same optimum found by exhaustive search.

The saved notebook output shows no reported discrepancies while runtimes for exhaustive search increase sharply with category cardinality. In the binary-target experiment, exhaustive search slows from well under a second at low cardinality to roughly 82 seconds total for `k=12` across 100 runs, while ordered splits remain much cheaper because they scale with the sorted category order rather than all partitions.

This is the key takeaway of the article and notebook: for common tree criteria such as Gini impurity in binary classification and MSE in regression, ordering categories by their target statistic gives a much more efficient search procedure than brute-force partitioning.

## How to run it

The notebook is self-contained and does not require any external dataset.

```bash
pip install pandas numpy matplotlib tqdm morethemes jupyter
cd mt_encoding
jupyter notebook mt_encoding.ipynb
```

## Notes

- The notebook uses synthetic data, so exact timings vary by machine.
- `morethemes` is only used for plot styling; the core experiment logic depends on `pandas`, `numpy`, `matplotlib`, and `tqdm`.
- The helper `mean_target_decode` function at the end is not used in the main experiment, but shows how an encoded mapping could be reversed when values are unique.
