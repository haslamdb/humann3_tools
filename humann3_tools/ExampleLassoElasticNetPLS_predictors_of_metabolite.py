import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.linear_model import Lasso, ElasticNet
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score

# ------------------------------------------------------------------
# 1. Generate Simulated Microbiome and Metabolite Data
# ------------------------------------------------------------------
np.random.seed(42)

n_samples = 200   # e.g., 200 mice
n_species = 1000  # e.g., 1000 microbial features

# Simulate microbiome data (X) - typically compositional
# Here, we'll just create random data for illustration.
X = np.random.gamma(shape=2.0, scale=1.0, size=(n_samples, n_species))

# OPTIONAL: Log-transform or CLR-transform if real compositional data
# For demonstration, let's do a simple log1p transform
X_log = np.log1p(X)

# Create a "true" relationship between some subset of features and y
true_features = 10
coeffs = np.zeros(n_species)
coeffs[:true_features] = np.linspace(5, 1, true_features)  # just an artificial pattern

# Generate the metabolite levels with some noise
y = X_log @ coeffs + np.random.normal(loc=0, scale=5, size=n_samples)

# ------------------------------------------------------------------
# 2. Split into Train and Test Sets
# ------------------------------------------------------------------
X_train, X_test, y_train, y_test = train_test_split(X_log, y, 
                                                    test_size=0.2, 
                                                    random_state=42)

# Optional: Standardize features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled  = scaler.transform(X_test)

# ------------------------------------------------------------------
# 3. LASSO Example
# ------------------------------------------------------------------
lasso = Lasso(random_state=42)

# Define hyperparameter grid for alpha (regularization strength)
param_grid_lasso = {
    'alpha': [0.01, 0.1, 1.0, 10.0, 100.0]
}

grid_lasso = GridSearchCV(
    estimator=lasso, 
    param_grid=param_grid_lasso,
    scoring='r2',        # use R^2 for regression
    cv=5,                # 5-fold cross-validation
    verbose=0,           # change to 1 or 2 for more logging
    n_jobs=-1            # use all available cores
)

grid_lasso.fit(X_train_scaled, y_train)
best_lasso = grid_lasso.best_estimator_

y_pred_lasso = best_lasso.predict(X_test_scaled)
r2_lasso = r2_score(y_test, y_pred_lasso)

print("=== LASSO Results ===")
print(f"Best alpha: {grid_lasso.best_params_['alpha']}")
print(f"Test R^2 : {r2_lasso:.3f}")
print(f"Number of non-zero coefficients: "
      f"{np.sum(best_lasso.coef_ != 0)}")
print()

# ------------------------------------------------------------------
# 4. Elastic Net Example
# ------------------------------------------------------------------
elastic = ElasticNet(random_state=42)

# Define hyperparameter grid for alpha and l1_ratio
# l1_ratio=1 is equivalent to LASSO, while l1_ratio=0 is Ridge
param_grid_enet = {
    'alpha': [0.01, 0.1, 1.0, 10.0],
    'l1_ratio': [0.2, 0.5, 0.8, 1.0]
}

grid_enet = GridSearchCV(
    estimator=elastic, 
    param_grid=param_grid_enet,
    scoring='r2',
    cv=5,
    verbose=0,
    n_jobs=-1
)

grid_enet.fit(X_train_scaled, y_train)
best_enet = grid_enet.best_estimator_

y_pred_enet = best_enet.predict(X_test_scaled)
r2_enet = r2_score(y_test, y_pred_enet)

print("=== Elastic Net Results ===")
print(f"Best parameters: {grid_enet.best_params_}")
print(f"Test R^2       : {r2_enet:.3f}")
print(f"Non-zero coefs : {np.sum(best_enet.coef_ != 0)}")
print()

# ------------------------------------------------------------------
# 5. Partial Least Squares (PLSRegression) Example
# ------------------------------------------------------------------
# PLS components must be <= min(n_samples, n_features). 
# Here, we guess a small number of components, but you can tune it.

# We'll do a simple grid search over number of components:
param_grid_pls = {
    'n_components': [2, 5, 10, 20]
}

pls = PLSRegression(scale=False)  # We already scaled X

grid_pls = GridSearchCV(
    estimator=pls,
    param_grid=param_grid_pls,
    scoring='r2',
    cv=5,
    verbose=0,
    n_jobs=-1
)

grid_pls.fit(X_train_scaled, y_train)
best_pls = grid_pls.best_estimator_

y_pred_pls = best_pls.predict(X_test_scaled)
r2_pls = r2_score(y_test, y_pred_pls)

print("=== PLS Regression Results ===")
print(f"Best n_components: {grid_pls.best_params_['n_components']}")
print(f"Test R^2         : {r2_pls:.3f}")

# Accessing the loadings, weights, or variable importance in PLS:
# - best_pls.x_weights_ : Weights for each X variable 
# - best_pls.x_loadings_: Loadings for each X variable
# - For feature importance, one can examine the absolute values of x_weights_ 
#   or use built-in methods from specialized PLS packages.

