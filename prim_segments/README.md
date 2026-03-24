# PRIM Segments

This subproject accompanies the article ["Find Unusual Segments in Your Data with Subgroup Discovery"](https://medium.com/data-science/find-unusual-segments-in-your-data-with-subgroup-discovery-2661a586e60c), also mirrored on [Towards Data Science](https://towardsdatascience.com/find-unusual-segments-in-your-data-with-subgroup-discovery-2661a586e60c/).

It contains a small, readable implementation of PRIM (Patient Rule Induction Method) for subgroup discovery and a notebook that applies it to a bank churn dataset.

## What is in this folder

- `prim.py`: a lightweight `PRIM` class for discovering high-response subgroups by iteratively peeling away low-quality regions of the feature space.
- `bank_churn_segments.ipynb`: a worked example on the "Churn for Bank Customers" dataset from Kaggle.

## What the code does

The implementation searches for segments where a binary target is unusually concentrated relative to a baseline rate. In this project, the target is customer churn (`exited`), and the algorithm looks for customer subgroups with elevated churn.

The `PRIM` class supports:

- numeric features, handled by quantile-based peeling on lower and upper bounds
- categorical features, handled by excluding category values one at a time
- repeated peeling steps that produce a sequence of increasingly specific boxes
- optional `patience`, which favors slightly larger but higher-value segments
- optional `bslevel`, which lets later iterations compare against a custom baseline

The segment score is based on the uplift of the subgroup target rate over the baseline, multiplied by subgroup size and optionally adjusted by `patience`.

## Reproduced example

The notebook uses the Kaggle dataset linked inside the notebook itself:

- dataset: "Churn for Bank Customers"
- target: `exited`
- feature set:
  `gender`, `geography`, `num_of_products`, `has_cr_card`, `is_active_member`,
  `credit_score`, `age`, `balance`, `tenure`, `estimated_salary`

The saved notebook outputs show three representative segments:

1. Default PRIM:
   `num_of_products < 2, is_active_member < 1, age > 37`
2. PRIM with `patience=2`:
   `num_of_products < 2, age > 37, age < 64`
3. A second-pass search after removing the first large churn segment:
   `gender != 'Male', num_of_products > 2, balance > 0.0`

These results match the main story of the article: adding patience can surface more useful segments than a purely greedy peeling strategy.

## How to run it

The notebook expects a `churn.csv` file in this directory.

Typical local setup:

```bash
pip install pandas numpy matplotlib jupyter
cd prim_segments
jupyter notebook bank_churn_segments.ipynb
```

If you want to use the class directly:

```python
import pandas as pd
from prim import PRIM

df = pd.read_csv("churn.csv")
X = df[[
    "gender", "geography", "num_of_products", "has_cr_card",
    "is_active_member", "credit_score", "age",
    "balance", "tenure", "estimated_salary"
]]
y = df["exited"]

model = PRIM(alpha=0.05, patience=2)
model.fit(X, y)

best = max(range(len(model.get_qual())), key=model.get_qual().__getitem__)
print(model.get_rules(best))
print(model.get_qual()[best])
print(model.get_size()[best])
```

## Notes and limitations

- `X` must be a pandas `DataFrame`.
- The target `y` is expected to be binary, with positive cases encoded as `1`.
- The implementation is intentionally compact and notebook-oriented rather than packaged as a general-purpose library.
- The notebook demonstrates one concrete dataset; it is easy to adapt the same logic to fraud, retention, quality, or risk segmentation problems.
