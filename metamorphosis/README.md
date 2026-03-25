# ML Metamorphosis

This subproject accompanies the article ["ML Metamorphosis: Chaining ML Models for Optimized Results"](https://medium.com/data-science/ml-metamorphosis-chaining-ml-models-for-optimized-results-d89d952627a9), also mirrored on [Towards Data Science](https://towardsdatascience.com/ml-metamorphosis-chaining-ml-models-for-optimized-results-d89d952627a9/).

It contains a small MNIST experiment whose main purpose is to improve a decision tree by first passing through a CNN stage. The CNN is not the final model of interest here; it is the intermediate model used to generate additional pseudo-labeled training data for the tree.

## What is in this folder

- `scripts/run_experiment.py`: entry point that runs the experiment, saves a timestamped result JSON under `output/`, and plots test accuracy by iteration.
- `self_training/experiment.py`: the main experiment loop, including MNIST loading, subset initialization, CNN self-training, pseudo-label selection, and train-set updates.
- `self_training/cnn.py`: the PyTorch CNN and its training helper.
- `self_training/experiment_result.py`, `self_training/metric_collection.py`, `self_training/plot.py`: result serialization, metric containers, and plotting helpers.
- `requirements.txt`: Python dependencies for the experiment.

## Main idea

The article's central idea is not just self-training in isolation. It is model chaining:

- Model A: train a CNN on a small labeled subset of MNIST
- Data 2: use the CNN's high-confidence predictions to create a larger pseudo-labeled transfer set
- Model B: train the target model, here a multiclass decision tree, on that enriched dataset

The intent is to get a better decision tree than you would get by fitting the tree directly on the limited labeled data alone. In the article, this is the key metamorphosis step: a stronger but less interpretable model helps produce data that a simpler, more interpretable model can learn from more effectively.

In this repository, the process is implemented on MNIST with:

- an initial labeled subset of `1000` training examples
- a confidence threshold of `0.98` for accepting pseudo-labels
- up to `10` self-training iterations

## What Was Added

This code is based on [Niklas von Moers' self-training project](https://github.com/NiklasvonM/Self-Training), but the main addition is the decision-tree stage that turns the self-training setup into the metamorphosis pipeline described in the article.

Specifically, this repository adds:

- training a multiclass `DecisionTreeClassifier` on the current training subset as it grows through pseudo-labeling
- evaluating that decision tree on the MNIST test set during each iteration
- the extra `scikit-learn` and `numpy` dependencies needed for that target-model step

## How to run it

From this directory:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_experiment.py
```

The script will:

- download MNIST into `./data` if needed
- run the iterative CNN self-training loop that expands the training subset
- retrain and evaluate the decision tree as that subset evolves
- save metrics to `output/experiment_result_<timestamp>.json`
- open a matplotlib plot of CNN test accuracy over iterations

## Notes

- The repository uses `requirements.txt`; the older `poetry` command that was previously in this README did not match the project files.
