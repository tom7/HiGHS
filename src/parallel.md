HiGHS currently has limited opportunities for exploiting parallel computing.

### Generally

By default, HiGHS runs in serial.

The parallelism set out below is
activated if the [parallel](https://ergo-code.github.io/HiGHS/options/definitions.html#parallel) option is set "on", by specifying
`--parallel` when running the
[executable](https://ergo-code.github.io/HiGHS/executable.html) via
the command line, or by setting it via a library call in an
application.

By default, when running in parallel, HiGHS will use half the
available threads on a machine. This number can be modified by setting
the value of the
[threads](https://ergo-code.github.io/HiGHS/options/definitions.html#threads)
option.

### Dual simplex

The dual simplex solver has a variant allowing concurrent
processing. This variant is used when the
[parallel](https://ergo-code.github.io/HiGHS/options/definitions.html#parallel)
option is set "on".

The concurrency used will be the value of
[simplex\_max\_concurrency](https://ergo-code.github.io/HiGHS/options/definitions.html#simplex_max_concurrency). If
this is fewer than the number of threads available, parallel
performance may be less than anticipated.

The speed-up achieved using the dual simplex solver is normally
bounded by the number of memory channels in the architecture, and
typically less than the values achieved by [Huangfu and
Hall](https://link.springer.com/article/10.1007/s12532-017-0130-5). This
is because enhancements to the serial dual simplex solver in recent
years have not been propagated to the parallel solver.

Unless an LP has significantly more variables than constraints, the
parallel dual simplex solver is unlikely to be worth using.

### MIP

The only parallel computation currently implemented in the MIP solver
occurs when performing symmetry analysis on the problem.

### Future plans

A prototype parallel LP solver has been developed, in which the
(serial) interior point solver and simplex variants are run
concurrently. When one runs to completion, the others are
stopped. However, to ensure that it runs deterministically requires
considerable further work. The non-deterministic solver will be
available by the end of 2023, but a deterministic solver is unlikely
to be available before the end of 2024.

The MIP solver has been written with parallel tree seach in mind, and
it is hoped that this will be implemented before the end of 2024. The
parallel LP solver will also enhance the MIP solver performance by
spoeeding up the solution of the root node.

Development of a parallel interior point solver will start in 2023,
and is expected to be completed by the end of 2024.



