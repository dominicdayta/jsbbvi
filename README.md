# Variance Control in Black Box Variational Inference Using the James-Stein Estimator

The following repository contains code and replication materials for the paper *Variance Control in Black Box Variational Inference Using the James-Stein Estimator* intended for submission to the 2024 Neural Information Processing System (NeurIPS) conference. We will continue to update this as the paper undergoes review and submission into other conferences/journals.

## Table of Contents

1. [Files Overview](#files-overview)
2. [Contributing](#contributing)

## Files Overview

The files are structured as follows:

1. `bbvi_simplemix.R` contains the base code and functions for performing BBVI (JS, Rao-Blackwellized, and Naive) on the simple gaussian mixture models, for use in the numerical experiments.
2. `variance_comparisons.ipynb` contains the numerical experiments demonstrating the variance reduction property of using the James-Stein estimator in BBVI.
3. `trace_plots.ipynb` contains the numerical experiments demonstrating that the James-Stein estimator converges to an optimal solution similar to (though at times not exactly) the Rao-Blackwellized solution.
4. `applied_model.ipynb` contains our application to a real dataset.

## Contributing

At the moment contributing to the code is not accepted.
