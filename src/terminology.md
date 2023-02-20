# HiGHS - High Performance Optimization Software
[![Build Status](https://github.com/ERGO-Code/HiGHS/workflows/build/badge.svg)](https://github.com/ERGO-Code/HiGHS/actions?query=workflow%3Abuild+branch%3Amaster)

Any linear optimization problem will have __decision variables__, a
linear or quadratic __objective function__, and linear __constraints__
and __bounds__ on the values of the decision variables. A
__mixed-integer__ optimization problem will require some or all of the
decision variables to take integer values. The problem may require the
objective function to be maximized or minimized whilst satisfying the
constraints and bounds. By default, HiGHS minimizes the objective
function.

### Bounds and the objective function

The bounds on a decision variable are the least and greatest values
that it may take, and infinite bounds can be specified. A linear
objective function is given by a set of coefficients, one for each
decision variable, and its value is the sum of products of
coefficients and values of decision variables. The objective
coefficients are often referred to as __costs__, and some may be
zero. When a problem has been solved, the optimal values of the
decision variables are referred to as the __(primal) solution__.

### Constraints and the feasible region

Linear constraints require linear functions of decision variables to
lie between bounds, and infinite bounds can be specified. If the
bounds are equal, then the constraint is an __equation__. If the
bounds are both finite, then the constraint is said to be __boxed__ or
__two-sided__. The set of points satisfying linear constraints and
bounds is known as the __feasible region__. Geometrically, this is a
multi-dimensional convex polyhedron, whose extreme points are referred
to as __vertices__.

### The constraint matrix

The coefficients of the linear constraints are naturally viewed as
rows of a __matrix__. The constraint coefficients associated with a
particular decision variable form a column of the constraint
matrix. Hence constraints are sometimes referred to as __rows__, and
decision variables as __columns__. Constraint matrix coefficients may
be zero. Indeed, for large practical problems it is typical for most
of the coefficients to be zero. When this property can be exploited to
computational advantage, the matrix is said to be __sparse__. When the
constraint matrix is not sparse, the solution of large problems is
normally intractable computationally.

### Optimization outcomes

It is possible to define a set of constraints and bounds that cannot
be satisfied, in which case the problem is said to be
__infeasible__. Conversely, it is possible that the value of the
objective function can be improved without bound whilst satisfying the
constraints and bounds, in which case the problem is said to be
__unbounded__. If a problem is neither infeasible, nor unbounded, it
has an __optimal solution__. The optimal objective function value for
a linear optimization problem may be achieved at more than point, in
which case the optimal solution is said to be __non-unique__.

### Basic solution of LP problems

An LP problem that is neither infeasible, nor unbounded, has an
optimal solution at a vertex. At a vertex, the decision variables can
be partitioned into as many __basic variables__ as there are
constraints, and __nonbasic variables__. Such a solution is known as a
__basic solution__, and the partition referred to as a __basis__.

### Dual values for continuous optimization

When none of the decision variables is required to take integer
values, the problem is said to be __continuous__. For
continuous problems, each variable and constraint has an
associated __dual variable__. The values of the dual
variables constitute the __dual solution__, and it is for
this reason that the term __primal solution__ is used to
distinguish the optimal values of the decision variables. At the
optimal solution of a continuous problem, some of the decision
variables and values of constraint functions will be equal to their
lower or upper bounds. Such a bound is said to
be __active__. If a variable or constraint is at a bound,
its corresponding dual solution value will generally be non-zero: when
at a lower bound the dual value will be non-negative; when at an upper
bound the dual value will be non-positive. When maximizing the
objective the required signs of the dual values are reversed. Due to
their economic interpretation, the dual values associated with
constraints are often referred to as __shadow prices__
or __fair prices__. Mathematically, the dual values
associated with variables are often referred to as __reduced
costs__, and the dual values associated with constraints are
often referred to as __Lagrange multipliers__.

### Sensitivity information for continuous optimization

Analysis of the change in optimal objective value of a continuous
linear optimization problem as the cost coefficients and bounds are
changed is referred to in HiGHS as __ranging__. For an
active bound, the corresponding dual value gives the change in the
objective if that bound is increased or decreased. This level of
analysis is often referred to as __sensitivity__. In
general, the change in the objective is only known for a limited range
of values for the active bound. HiGHS will return the limits of
these __bound ranges__ ranges, the objective value at
both limits and the index of a variable or constraint that will
acquire an active bound at both limits. For each variable with an
active bound, the solution will remain optimal for a range of values
of its cost coefficient. HiGHS will return the values of
these __cost ranges__. For a variable or constraint whose
value is not at a bound, HiGHS will return the range of values that
the variable or constraint can take, the objective values at the
limits of the range, and the index of a variable or constraint with a
bound that will become in active at both limits.



