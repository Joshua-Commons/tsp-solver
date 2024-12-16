# Metric TSP Solver

This codebase contains a heuristic to find solutions to metric TSP problems which are of size 10-50, within 30 seconds.

## Table of Contents

- [Metric TSP Solver](#metric-tsp-solver)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
  - [Requirements](#requirements)
  - [Theory](#theory)
    - [Held-Karp Algorithm](#held-karp-algorithm)
    - [Christofides Algorithm](#christofides-algorithm)
  - [Usage](#usage)

## Installation

Copy the GitHub repo:
```
    git clone https://github.com/Joshua-Commons/tsp-solver.git
```

## Requirements

Python version `3.10.16` and NumPy version `2.1.3`. 

Optionally, load conda environment `tsp_env` by `conda env create -f tsp_env.yml`.

Activate environment with `conda activate tsp_env`.

## Theory

### Held-Karp Algorithm

The Held-Karp algorithm uses dynamic programming with bit masking to solve a TSP. This algorithm is not a heuristic, and so returns the exact solution. [Wikipedia Link](!https://en.wikipedia.org/wiki/Held%E2%80%93Karp_algorithm).

The time complexity of the algorithm is `O(n^2 * 2^n)` and is so for our purposes is only feasible for graphs which are up to size `n=23`.

### Christofides Algorithm

The Christofides Algorithm generates a minimum spanning tree and performs modifications to it to make a valid TSP solution. [Wikipedia Link](!https://en.wikipedia.org/wiki/Christofides_algorithm).

## Usage

1. Initialise object with `solver = TSP_Solver()`

2. Load a distance matrix by either:

    1. Generating a distance matrix by `matrix = solver.generate_matrix(n)` where 
    `n` is the desired size of matrix 10-50 inclusive.

    2. Load a distance matrix from a `.csv` file by `matrix = solver.load_matrix(path)` where
    `path` is the path to your file. Ideally, under the `data` folder.

3. Validate matrix form (done automatically when `solver.solve(matrix)` is called.)

4. Find solution to metric TSP with `solver.solve(matrix)` where `matrix` is the
distance matrix.