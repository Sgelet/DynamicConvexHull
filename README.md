# Simple and Robust Dynamic Two-Dimensional Convex Hull

This is the repository and testbed used for experimental evaluation as presented in the paper "Simple and Robust Dynamic Convex Hull".
It contains the header required to use the data structures, scripts to generate synthetic data as well as perform evaluations of performance in dynamic scenarios. Included is also a rudimentary tool to replicate plots as shown in the paper.

## Run options

To manually run tests, simply build and run with the following program arguments, and supply the data through standard in:

- `RUN X` runs the dynamic scenario, reporting every X operations. Reads input from stdin.
- `GEN X Y` generates X numbers according to the distribution specified by Y as an index into the following list of distributions: \[normal, uniform, uniform cropped to circle, circle\]. Outputs the generated points to stdout.
- `VER X` runs verification against CGAL static hulls, comparing every X updates. Reads input from stdin.

The script `smalltest.sh` shows an example of how to automate the process.

## Plots

We provide a very rudimentary tool for generating plots through the Python program `output_parser.py`.