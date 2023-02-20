### Specification

HiGHS is software for definition, modification and solution of general bounded large scale sparse linear optimization problems. It can solve linear programming (LP) problems of the form
```math
\textrm{min} \qquad c^Tx \qquad \textrm{subject to} \qquad L \le Ax \le U; \qquad l \le x \le u.
```
and mixed integer programming (MIP) problems of the same form, for
which some of the variables must take integer values. HiGHS also
solves quadratic programming (QP) problems, with objective term
$(1/2)x^TQx$, where the Hessian matrix $Q$ is positive
semi-definite. It cannot solve QP problems where some of the variables
must take integer values.

More on the
[terminology](http://ergo-code.github.io/HiGHS/terminology.html) of
optimization is available.

### Solvers

For LP problems, HiGHS has implementations of both the revised simplex
and interior point methods. MIPs are solved by branch-and-price, and
QPs by active set.

### Model and solution management

HiGHS has comprehensive tools for defining and extracting models. This
can be done either to/from MPS or (CPLEX) format LP files, or via
method calls. HiGHS also has methods that permit the incumbent model
to be modified. Soluitons can be supplied and extracted using either
files or method calls.

### OS

HiGHS can be used on Windows, Linux and MacOS.

### Compilers

HiGHS can be used with the following compilers:

- Clang ` clang `
- GNU ` g++ ` 
- Intel ` icc `

### Dependencies

No third party sortware is required by HiGHS, except for the Threads library.

In order to build HiGHS from source CMake 3.15 is required. For precompiled executables and libraries please contact us at [highsopt@gmail.com](mailto:highsopt@gmail.com).

### Reference

[Parallelizing the dual revised simplex method](http://www.maths.ed.ac.uk/hall/HuHa13/)
Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: 10.1007/s12532-017-0130-5

#### Performance

The performance of HiGHS relative to some commercial and open-source simplex solvers may be assessed via the [Mittelmann benchmarks](http://plato.asu.edu/ftp/lpsimp.html).
