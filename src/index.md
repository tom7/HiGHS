# HiGHS - High Performance Optimization Software
[![Build Status](https://github.com/ERGO-Code/HiGHS/workflows/build/badge.svg)](https://github.com/ERGO-Code/HiGHS/actions?query=workflow%3Abuild+branch%3Amaster)

HiGHS is software for the definition, modification and solution of
large scale sparse linear optimization problems. It is freely
available from [GitHub](https://github.com/ERGO-Code/HiGHS) under the
MIT licence, and has no third-party dependencies.

HiGHS can be used as a standalone
[executable](https://ergo-code.github.io/HiGHS/executable.html) on
Windows, Linux and MacOS. There is also a C++11 library that can be
used within a C++ project or, via its C, C#, FORTRAN, Julia and Python
interfaces.

[Building HiGHS from source
code](https://ergo-code.github.io/HiGHS/cpp/get-started.html#Building-HiGHS-from-source-code)
requires a C++ compiler and CMake, but no other third-party utilities.

Precompiled binaries are available for a variety of platforms at the
[JuliaBinaryWrappers HiGHS
repository](https://github.com/JuliaBinaryWrappers/HiGHS_jll.jl/releases). Each
includes library files for linking to external projects, and a
stand-alone executable

For Windows users: if in doubt, choose the `x86_64-w64-mingw32-cxx11.tar.gz` file

For Mac users: choose the `x86_64-apple-darwin.tar.gz` file.

### Specification

HiGHS can solve linear programming (LP) problems of the form
```math
\textrm{minimize} \qquad c^Tx \qquad \textrm{subject to} \qquad L \le Ax \le U; \qquad l \le x \le u.
```
and mixed integer programming (MIP) problems of the same form, for
which some of the variables must take integer values. HiGHS also
solves quadratic programming (QP) problems, with (additional)
objective term $(1/2)x^TQx$, where the Hessian matrix $Q$ is positive
semi-definite. It cannot solve QP problems where some of the variables
must take integer values.

More on the
[terminology](http://ergo-code.github.io/HiGHS/terminology.html) of
optimization is available.

### Solvers

For LP problems, HiGHS has implementations of both the revised simplex
and interior point methods. MIPs are solved by branch-and-price, and
QPs by active set.

###  Reference

Academic users of HiGHS should cite the following article.

[Parallelizing the dual revised simplex method](http://www.maths.ed.ac.uk/hall/HuHa13/)
Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: 10.1007/s12532-017-0130-5

### Performance

The performance of HiGHS relative to some commercial and open-source simplex solvers may be assessed via the [Mittelmann benchmarks](http://plato.asu.edu/ftp/lpsimp.html).


### Feedback

Your comments or specific questions on HiGHS would be greatly
appreciated, so please send an email to
[highsopt@gmail.com](mailto:highsopt@gmail.com) to get in touch with
the team.
