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

Here, we present some of the main results from the analyses performed with this code package.  We will use three different efficiency functions.

    1. A longitude independent efficiency function, supplied directly in [1705.00009](https://arxiv.org/pdf/1705.00009.pdf).  We will label this as `long-indep`
    2. A longitude-corrected efficiency function, corrected using data in [1305.4385](https://arxiv.org/abs/1305.4385).  We will label this as (`2PC-corrected`)
    3. A longitude-corrected efficiency function, corrected by rescaling the longitude-independent efficiency function by the actual data values within the pixels.  We will label this as (`data-corrected`)

The reason we use three different efficiency functions is that this allows us to understand the systematic uncertainty associated with the fact that we likely do not have the correct efficiency.

Below, unless floated, values are taken to be:

- $n = 2.35$,
- $\sigma = 1.528$,
- $\alpha = 2.6$,
- $\beta = 1.2$,
- $z_0 = 0.7$.

Priors used: 

- $N_B, N_D$: [0,3000000]
- $\alpha$: [2.1, 5.0]
- $z_0$: [0.01, 2.0]
- $\beta$: [1.1,3.0]

Floating parameters as in Table 2 of [[1705.00009](https://arxiv.org/pdf/1705.00009.pdf):



All results below are computed using `multinest`.  



<!-- ### `Long-Indep`

| $N_D$             | $N_B$                     | $z_0$         | $\beta$       | $\alpha$      | TS |
|-------------------|---------------------------|---------------|---------------|---------------|----|
| $1483854^{+555529}_{-453093}$ | $0$                         | $0.08^{+0.05}_{-0.03}$ | $2.09^{+0.07}_{-0.08}$    | -             | $0$  |
|$1153102^{+525377}_{-378358}$      | $506568^{+496453}_{-238196}$    |  $0.05^{+0.03}_{-0.02}$| $2.10^{+0.08}_{-0.07}$    | $2.6$             | $10.34$ |
| $1052708^{+456557}_{-343513}$ | $966064^{+1160994}_{-587532}$ | $0.06^{+0.04}_{-0.02}$ | $2.08^{+0.07}_{-0.07}$ | $2.83^{+0.11}_{-0.27}$ | $11.39$ |


#### Float $N_D$, $z_0$ and $\beta$

![Disk only](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-z0-beta-lon-indep_eff.png "Disk only")

#### Float $N_D$, $N_B$, $z_0$ and $\beta$

![Disk and bulge](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-nb-z0-beta-lon-indep_eff.png "Disk and bulge")

#### Float $N_D$, $N_B$, $z_0$, $\beta$ and $\alpha$

![Float alpha](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-nb-alpha-z0-beta-lon-indep_eff.png "Float alpha")


### `2PC-Corrected`

| $N_D$             | $N_B$                     | $z_0$         | $\beta$       | $\alpha$      | TS |
|-------------------|---------------------------|---------------|---------------|---------------|----|
| $1482073^{+575204}_{-465966}$ | $0$                         | $0.08^{+0.05}_{-0.03}$ | $2.09_{-0.08}^{+0.08}$    | -             | 0  |
|$1068731^{+483159}_{-355819}$      | $590663^{+540377}_{-282268}$    |  $0.05^{+0.04}_{-0.02}$| $2.10^{+0.08}_{-0.08}$    | $2.6$             | $13.21$ |
| $955272^{+409682}_{-318262}$ | $1328045.71^{+1066159}_{-803946}$ | $0.06^{+0.04}_{-0.02}$ | $2.08^{+0.07}_{-0.07}$ | $2.86^{+0.08}_{-0.20}$ | $15.28$ |


#### Float $N_D$, $z_0$ and $\beta$

![Disk only](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-z0-beta.png "Disk only")

#### Float $N_D$, $N_B$, $z_0$ and $\beta$

![Disk and bulge](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-nb-z0-beta.png "Disk and bulge")

#### Float $N_D$, $N_B$, $z_0$, $\beta$ and $\alpha$

![Float alpha](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-nb-alpha-z0-beta.png "Float alpha")

 
### `data-corrected`

| $N_D$             | $N_B$                     | $z_0$         | $\beta$       | $\alpha$      | TS |
|-------------------|---------------------------|---------------|---------------|---------------|----|
| $1516977^{+559102}_{-445328}$ | $0$                         | $0.07^{+0.05}_{-0.03}$ | $2.09^{+0.07}_{-0.07}$    | -             | $0$  |
|$1117715^{+484416}_{-363228}$      | $533365^{+529982}_{-246334}$    |  $0.05^{+0.03}_{-0.02}$| $2.10^{+0.08}_{-0.07}$    | $2.6$             | $11.19$ |
| $1022473^{+463709}_{-331367}$ | $1120457^{+1080567}_{-671964}$ | $0.06^{+0.04}_{-0.02}$ | $2.08^{+0.08}_{-0.07}$ | $2.84^{+0.09}_{-0.25}$ | $12.84$ |


#### Float $N_D$, $z_0$ and $\beta$

![Disk only](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-z0-beta-data_eff.png "Disk only")

#### Float $N_D$, $N_B$, $z_0$ and $\beta$

![Disk and bulge](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-nb-z0-beta-data_eff.png "Disk and bulge")

#### Float $N_D$, $N_B$, $z_0$, $\beta$ and $\alpha$

![Float alpha](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-nb-alpha-z0-beta-data_eff.png "Float alpha")

  -->
