# GCE-2FIG

[![arXiv](https://img.shields.io/badge/arXiv-1710.10266%20-green.svg)](https://arxiv.org/abs/1710.10266)

This code implements the analysis described in [1710.10266](https://arxiv.org/pdf/1710.10266.pdf) and uses 2FIG sources to constrain the spatial distributions of pulsar-like PS populations. The main purpose of the code is to test the results presented in [1705.00009](https://arxiv.org/pdf/1705.00009.pdf).

All of the code is contained in the `likelihood` directory and analysis notebooks are contained in `notebooks/`.

## Authors

-  Richard Bartels; richard.t.bartels at gmail dot com
-  Dan Hooper; dhooper at fnal dot gov
-  Tim Linden; linden.70 at osu dot edu
-  Siddharth Mishra-Sharma; smsharma at princeton dot edu
-  Nicholas Rodd; nrodd at mit dot edu
-  Benjamin Safdi; bsafdi at umich dot edu
-  Tracy Slatyer; tslatyer at mit dot edu

## Compiling and running

This code package is written in `python` and `cython`. To compile the `cython`, go into the `likelihood/` subfolder and execute the file `make.sh`.  Note that you must have `cython` installed.  
Examples of how to run the code are presented in the `run/` folder.

## Parameter scans

The scripts to run the scans in parallel on the cluster are located in the `run/` folder, also containing examples of `SLURM` batch files used to launch these scan. The file names are below. 
These are used to obtain the Monte Carlo chains, which are included in the repository in `run/chains/`.

Three scans are performed, floating:

- ND, z0 and beta (`run_nd_z0_beta_prior.py`)
- ND, NB, z0 and beta (`run_nd_nb_z0_beta_prior.py`)
- ND, NB, z0, beta and alpha (`run_nd_nb_z0_beta_alpha_prior.py`).

The likelihoods of the runs including a bulge population are compared with one without it to test the preference for a bulge population of pulsars. 
This is analyzed in the notebook `Full_Scan.ipynb` as described below.

## Notebooks

The following notebooks are used in the analysis:

1. `Efficiency.ipynb`: Constructing of the pulsar detection efficiency for binned data.
2. `Sources.ipynb`: Constructing the data used in the analysis. The binning scheme is described and data is constructed from a combination of 2FIG candidates and 3FGL sources.
3. `Models_and_Likelihoods.ipynb`: Deriving the modeled number of disk and bulge sources in each bin. Calculating the likelihood for a given set of disk and bulge model parameters.
4. `Full_Scan.ipynb`: Analyzing the results of the full parameter scans above, extracting best-fit parameters and TS in preference of a bulge population.
5. `A_Resolved_Fraction.ipynb`: Calculating the fraction of the GCE that should be resolved for various luminosity functions, including the one calculated in the original 2FIG paper.

## Results

The best fit parameters taken from Table 2 of [1705.00009](https://arxiv.org/pdf/1705.00009.pdf) are given as:

![alt text](https://github.com/bsafdi/GCE-2FIG/blob/master/notebooks/plots/Table_Fermi.png "Fermi best fit parameters")

Our runs indicate instead the following best fit parameters:

![alt text](https://github.com/bsafdi/GCE-2FIG/blob/master/notebooks/plots/Table_Us.png "Our best fit parameters")

In particular our returned TS values, whilst still in favor of a bulge population, are more modest. Further we favor different values for both z0 and beta. The derivation of these parameters is shown in the notebooks outlined below.
