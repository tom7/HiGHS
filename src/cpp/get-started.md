## Download source code

HiGHS can be cloned from the Edinburgh Group in Research and Optimization ([ERGO](https://www.maths.ed.ac.uk/ERGO/)) [GitHub repo](https://www.github.com/ERGO-COde/HiGHS).

``` bash
git clone https://github.com/ERGO-Code/HiGHS.git
```


### Build HiGHS from source code

HiGHS uses CMake as a build system. The simplest setup is to create a build folder (within the folder into which HiGHS has been downloaded) and then build HiGHS within it. The name of the build folder is arbitrary but, assuming it is HiGHS/build, the full sequence of commands required is as follows

``` bash
cd HiGHS
mkdir build
cd build
cmake -DFAST_BUILD=ON ..
cmake --build . 
```

This creates the executable `build/bin/highs`.

### Test Build

To perform a quick test to see whether the compilation was successful, run ctest from within the build folder.

``` bash
ctest 
```

### Install 

The default installation location may need administrative permissions. To install, after building and testing, run 

``` bash
cmake --install . 
```

To install in a specified installation directory run CMake with the `CMAKE_INSTALL_PREFIX` flag set: 

``` bash
cmake -DFAST_BUILD=ON -DCMAKE_INSTALL_PREFIX=/path/to/highs_install ..
cmake --build .
cmake --install . 
```

## Precompiled executables 
Precompiled executables are available for a variety of platforms at the [JuliaBinaryWrappers HiGHS repository](https://github.com/JuliaBinaryWrappers/HiGHS_jll.jl/releases).

Note that HiGHS is still pre-1.0, so the version numbers in the releases do not match versions of HiGHS in this repository.

For Windows users: if in doubt, choose the `x86_64-w64-mingw32-cxx11.tar.gz` file

For Mac users: choose the `x86_64-apple-darwin.tar.gz` file.

### Input file formats

HiGHS can parse .mps and .lp files. Models can also be loaded at runtime from another program using the library interface.

## Running the executable

Assuming the executable was created following the [build]() instructions, and the lp model is specified in `model.mps` (see more on LP [Input file formats])In the following discussion, the name of the executable file created in build/bin when building HiGHS is assumed to be highs. HiGHS can read plain text MPS files and LP files (but not compressed files), and the following command solves the model in model.mps

```bash
  /build/bin/highs /path/to/model.mps
```

For usage and runtime option information see 

``` bash
./bin/highs --help
```
### Command line options

When HiGHS is run from the command line, some fundamental option values may be specified directly. Many more may be specified via a file. Formally, the usage is

```bash
./build/bin/highs --help
Running HiGHS 1.0.0
Copyright (c) 2021 ERGO-Code under MIT licence terms

HiGHS options
Usage:
  ./build/bin/highs [OPTION...] [file]

      --model_file arg    File of model to solve.
      --presolve arg      Presolve: "choose" by default - "on"/"off" are
                          alternatives.
      --solver arg        Solver: "choose" by default - "simplex"/"ipm" are
                          alternatives.
      --parallel arg      Parallel solve: "choose" by default - "on"/"off"
                          are alternatives.
      --time_limit arg    Run time limit (double).
      --options_file arg  File containing HiGHS options.
  -h, --help              Print help.
```