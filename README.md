# GCE-2FIG
This code uses 2FIG sources to constrain spatial distributions of pulsar-like PS populations.  The main purpose of the code is to test the results presented in [1705.00009](https://arxiv.org/pdf/1705.00009.pdf).

## Compiling and running

This code package is written in `python` and `cython`. To compile the `cython`, go into the `likelihood/` subfolder and execute the file `make.sh`.  Note that you must have `cython` installed.  Examples of how to run the code are presented in the `examples/` subfolder.

## Notebooks

The following notebooks are used in the analysis:

1. `Efficiency.ipynb`: Constructing of the pulsar detection efficiency for binned data in several ways.
2. `Sources.ipynb`: Constructing the data used in the analysis. The binning scheme is described and data is constructed from 2FIG candidates and 3FGL sources.
3. `Models and Likelihoods.ipynb`: Deriving the modeled number of disk and bulge sources in each bin. Calculating the likelihood for a given set of disk and bulge model parameters.
4. `Simple Scans.ipynb`: Performing simple scans with a limited number of parameters and reduced accuracy for illustration.
5. `Full Scans and Results.ipynb`: Performing the full scan over model parameters and analyzing the results.

## Results

The best fit parameters taken from Table 2 of [1705.00009](https://arxiv.org/pdf/1705.00009.pdf) are given as:

![alt text](https://github.com/bsafdi/GCE-2FIG/blob/master/notebooks/Table_Fermi.png "Fermi best fit parameters")

Our runs indicate instead the following best fit parameters:

![alt text](https://github.com/bsafdi/GCE-2FIG/blob/master/notebooks/Table_Fermi.png "Our best fit parameters")
