import numpy as np
import pandas as pd

class PRIM:
    def __init__(self, alpha=0.05, patience=0, bslevel=None):
        self.alpha = alpha  # peeling parameter
        self.patience = patience  # extra parameter to control patience
        self.bslevel = bslevel  # for churn data: benchmark level of churn_rate
        # (https://medium.com/towards-data-science/figuring-out-the-most-unusual-segments-in-data-af5fbeacb2b2)
        self.boxes_ = []  # To store boxes
        self.qualities_ = []  # To store qualities of boxes
        self.sizes_ = []  # To store relative sizes of boxes
        self.column_names_ = None  # To store column names of input data

    def _target_fun(self, npos, n):
        bsl = self.Np_ / self.N_ if self.bslevel is None else self.bslevel
        tar = n * (npos / n - bsl)
        return tar

    def _box_quality(self, in_box):
        # Count points in the trial box
        n_in_box = np.count_nonzero(in_box)
        # Calculate points labelled 1 in the trial box
        n_pos_in_box = np.count_nonzero(self.y_[in_box])
        # Calculate quality
        quality = self._target_fun(n_pos_in_box, n_in_box)*(n_in_box**self.patience)
        return quality, n_in_box

    def _upd_in_iteration_vars(self, in_box, n_in_box, quality):
        self.ind_in_best_box_ = in_box
        self.n_in_best_box_ = n_in_box
        self.max_quality_ = quality
        return self

    def _peel_num(self, idx, column):
        column_data = self.X_[column]
        for direction in [0, 1]:  # 0 for lower bound, 1 for upper bound
            if direction == 0:
                cut_val = column_data.quantile(self.alpha, interpolation='lower')
                if cut_val == max(column_data):
                    in_box = column_data >= cut_val
                else:
                    in_box = column_data > cut_val
            else:
                cut_val = column_data.quantile(1 - self.alpha, interpolation='higher')
                if cut_val == min(column_data):
                    in_box = column_data <= cut_val
                else:
                    in_box = column_data < cut_val

            quality, n_in_box = self._box_quality(in_box)
            if quality > self.max_quality_:
                self._upd_in_iteration_vars(in_box, n_in_box, quality)
                self.best_box_ = [b.copy() for b in self.box_]
                self.best_box_[idx][direction] = cut_val
        return self

    def _peel_cat(self, idx, column):
        column_data = self.X_[column]
        for value in column_data.unique():
            in_box = column_data != value
            quality, n_in_box = self._box_quality(in_box)
            if quality > self.max_quality_:
                self._upd_in_iteration_vars(in_box, n_in_box, quality)
                self.best_box_ = [b.copy() for b in self.box_]
                self.best_box_[idx].append(value)
        return self

    def _set_fit_vars(self):
        self.N_ = len(self.y_)
        self.Np_ = np.sum(self.y_)
        self.box_ = self._get_initial_restrictions(self.X_)
        self.column_names_ = self.X_.columns
        return self

    def _set_in_iteration_vars(self):
        self.max_quality_ = -np.inf
        self.best_box_ = None
        self.n_in_best_box_ = self.N_
        self.ind_in_best_box_ = None
        return self

    def _upd_fit_vars(self):
        # Update the box and record it
        self.box_ = self.best_box_
        self.boxes_.append(self.box_)
        self.qualities_.append(self.max_quality_/(self.n_in_best_box_**self.patience))
        self.sizes_.append(self.n_in_best_box_)
        self.X_ = self.X_[self.ind_in_best_box_]
        self.y_ = self.y_[self.ind_in_best_box_]
        return(self)

    def fit(self, X, y):
        # Ensure X is a DataFrame
        if not isinstance(X, pd.DataFrame):
            raise ValueError("X must be a pandas DataFrame")

        self.X_ = X
        self.y_ = y
        self._set_fit_vars()

        for iteration in range(120):
            self._set_in_iteration_vars()
            for idx, column in enumerate(self.X_.columns):
                column_data = self.X_[column]
                if column_data.nunique() > 1:
                    if column_data.dtype.kind in 'bifc':  # numeric column
                        self._peel_num(idx, column)
                    else:                                 # categorical column
                        self._peel_cat(idx, column)

            # If we can't find a suitable box or the box contains less than 1/alpha points, break
            if self.best_box_ is None or self.n_in_best_box_ < 1 / self.alpha:
                break

            # Update the box and record it
            self._upd_fit_vars()
        return self

    def _get_initial_restrictions(self, X):
        box = []
        for column in X.columns:
            if X[column].dtype.kind in 'bifc':  # Check if the column is numeric (boolean, integer, float, complex)
                # For numeric data, use range as initial restriction
                box.append([-np.inf, np.inf])
            else:
                # For non-numeric data, use unique values as initial restriction
                box.append([])
        return box

    def get_nrestr(self):
        return [np.count_nonzero(~np.isinf(box)) for box in self.boxes_]

    def get_rules(self, ind):
        rules = []
        for box, clnm in zip(self.boxes_[ind], self.column_names_):
            if not box:
                continue  # Skip empty lists
            elif isinstance(box[0], str):
                # Handling strings (including multiple strings)
                conditions = [f"{clnm} != '{value}'" for value in box]
                rules.append(' or '.join(conditions))
            else:
                # Handling numeric ranges
                if box[0] == float('-inf') and box[1] != float('inf'):
                    rules.append(f"{clnm} < {box[1]}")
                elif box[0] != float('-inf') and box[1] == float('inf'):
                    rules.append(f"{clnm} > {box[0]}")
                elif box[0] != float('-inf') and box[1] != float('inf'):
                    rules.append(f"{clnm} > {box[0]}")
                    rules.append(f"{clnm} < {box[1]}")
                # Add additional conditions if needed for other numeric range scenarios
        return ', '.join(rules)

    def get_qual(self):
        return self.qualities_

    def get_box(self, ind):
        return self.boxes_[ind]

    def get_size(self):
        return self.sizes_
