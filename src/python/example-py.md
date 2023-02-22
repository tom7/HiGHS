
If highspy is not already installed, run

```
pip install highspy
```

## Initialize Highs

```
import highspy
import numpy as np

# Highs h
h = highspy.Highs()
```

## Load a model

```
# Load a model from MPS file model.mps
filename = 'model.mps'
h.readModel(filename)
```

## Build a model

Build the model

Firstly, via a sequence of calls to __addVar__ and __addRow__.
```
# Define two variables, first using identifiers for the bound values,
# and then using constants
lower = 0
upper = 4
h.addVar(lower, upper)
h.addVar(1, inf)

# Define the objective coefficients (costs) of the two variables,
# identifying the variable by index, and changing its cost from the
# default value of zero
cost = 1
h.changeColCost(0, cost)
h.changeColCost(1, 1)

# Define constraints for the model
#
# The first constraint (x_1<=7) has only one nonzero coefficient,
# identified by variable index 1 and value 1
num_nz = 1
index = 1
value = 1
h.addRow(-inf, 7, num_nz, index, value)

# The second constraint (5 <= x_0 + 2x_1 <= 15) has two nonzero
# coefficients, so arrays of indices and values are required
num_nz = 2
index = np.array([0, 1])
value = np.array([1, 2])
h.addRow(5, 15, num_nz, index, value)

# The final constraint (6 <= 3x_0 + 2x_1) has the same indices but
different values num_nz = 2 value = np.array([3, 2]) h.addRow(6, inf,
num_nz, index, value)

# Access LP 
lp = h.getLp()
num_nz = h.getNumNz()
print('LP has ', lp.num_col_, ' columns', lp.num_row_, ' rows and ', num_nz, ' nonzeros')

```

## Pass a model

```
# Pass a model from a HighsLp instance
inf = highspy.kHighsInf
lp = highspy.HighsLp()
lp.num_col_ = 2
lp.num_row_ = 2
lp.sense_ = highspy.ObjSense.kMaximize
lp.col_cost_ = np.array([8, 10], dtype=np.double)
lp.col_lower_ = np.array([0, 0], dtype=np.double)
lp.col_upper_ = np.array([inf, inf], dtype=np.double)
lp.row_lower_ = np.array([-inf, -inf], dtype=np.double)
lp.row_upper_ = np.array([120, 210], dtype=np.double)
lp.a_matrix_.start_ = np.array([0, 2, 4])
lp.a_matrix_.index_ = np.array([0, 1, 0, 1])
lp.a_matrix_.value_ = np.array([0.3, 0.7, 0.5, 0.5], dtype=np.double)
h.passModel(lp)
```

## Solve problem

```
h.run()
```

## Print solution information 
```
solution = h.getSolution()
basis = h.getBasis()
info = h.getInfo()
model_status = h.getModelStatus()
print('Model status = ', h.modelStatusToString(model_status))
print()
print('Optimal objective = ', info.objective_function_value)
print('Iteration count = ', info.simplex_iteration_count)
print('Primal solution status = ', h.solutionStatusToString(info.primal_solution_status))
print('Dual solution status = ', h.solutionStatusToString(info.dual_solution_status))
print('Basis validity = ', h.basisValidityToString(info.basis_validity))
```