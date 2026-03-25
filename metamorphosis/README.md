# ML Metamorphosis

This subproject explores a simple "ML metamorphosis" workflow on MNIST: start with a small labeled subset, train a neural network, use the network to pseudo-label high-confidence unlabeled examples, add those examples back into the training set, and repeat. The idea matches the article referenced in [`AGENT_TASK.md`](/workspaces/medium/metamorphosis/AGENT_TASK.md): one model's output changes the effective training data for the next training round.

The codebase is adapted from Niklas von Moers' self-training experiment and extends it with an additional classical model step: training and evaluating a multiclass `DecisionTreeClassifier` on the same evolving training subset before each CNN self-training iteration.

## Project layout

- [`scripts/run_experiment.py`](/workspaces/medium/metamorphosis/scripts/run_experiment.py) is the entry point. It configures the experiment, runs several iterations, saves metrics to `output/*.json`, and plots test accuracy over time.
- [`self_training/experiment.py`](/workspaces/medium/metamorphosis/self_training/experiment.py) contains the core experiment loop: MNIST loading, labeled subset initialization, CNN training, pseudo-label generation, train-set expansion, and the added decision-tree training and evaluation.
- [`self_training/cnn.py`](/workspaces/medium/metamorphosis/self_training/cnn.py) defines the PyTorch CNN and its training helper.
- [`self_training/experiment_result.py`](/workspaces/medium/metamorphosis/self_training/experiment_result.py), [`self_training/metric_collection.py`](/workspaces/medium/metamorphosis/self_training/metric_collection.py), and [`self_training/plot.py`](/workspaces/medium/metamorphosis/self_training/plot.py) handle metric serialization and visualization.
- [`requirements.txt`](/workspaces/medium/metamorphosis/requirements.txt) lists the Python dependencies used by the experiment.

## Running the experiment

From [`metamorphosis/`](/workspaces/medium/metamorphosis):

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_experiment.py
```

What the run does:

- downloads the MNIST dataset into `./data` if it is not already present
- initializes the experiment with `initial_subset_size=1000` and `confidence_threshold=0.98`
- runs up to 10 self-training iterations
- writes a timestamped result JSON under `output/`
- opens a matplotlib plot of test accuracy by iteration

## What Was Added

Compared with the original self-training codebase, this version adds:

- decision-tree training on the current labeled subset via scikit-learn
- decision-tree evaluation on the MNIST test set each iteration
- the extra scikit-learn and NumPy dependencies needed for that classical-model baseline step

## Notes

- The current experiment loop prints decision-tree accuracy, but only the CNN self-training metrics are persisted to the saved experiment result.
- The repository currently uses plain `pip` plus [`requirements.txt`](/workspaces/medium/metamorphosis/requirements.txt); the old `poetry run ...` instruction in the previous README did not match the project files.
