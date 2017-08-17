# GCE-2FIG
This code uses 2FIG sources to constrain spatial distributions of pulsar-like PS populations.  The main purpose of the code is to test the results presented in 1705.00009.

## Compiling and Running

This code package is written in `python` and `cython`.  To compile the `cython`, go into the `python/` subfolder and execute the file `make.sh`.  Note that you must have `cython` installed.  Examples of how to run the code are presented in the `examples/` subfolder.

## The likelihood Function

The `GCE-2FIG` code package uses the likelihood function presented in 1705.00009.  The likelihood is imlemented in the file `likelihood.pyx` in the `python/` subfolder.  The sky is spatially binned in a cartesian grid, and the data consists of the number of pulsar candidates detected in each bin.  The model parameters characterize the disk and Bulge point source populations, and as we scan over the model parameters we calculated the expected number of detected sources in each bin, utlizing the _Fermi_-provided efficiency function for detecting sources at a given flux and spatial position.  The likelihood is then given by the product over all pixels of the Poisson probabilities to observe the data counts given the predicted model counts.

As an option, we also incorporate the prior distribution in 1705.00009.

## The Expected Number PS Counts

Here, we describe how we calculate the expected number of PS counts, given the model parameters.  These computations are implemented in the file `get_counts_inline.pyx` in the `python/` subfolder.  First, we describe the computation for a general PS population, with density function $\rho(s,\ell,b)$, where $(s,\ell,b)$ are spatial coordinates with origin at the Earth, and luminosity function $dN/dL$.  Note that $\ell$ and $b$ are Galactic longitude and latitude, respectively, while $s$ is the distance from the Earth.  The expected number of counts in a given spatial bin (labeled by $i,j$) and flux bin (labelled by $k$) is given by the integral

$$
N_{i,j,k} = \omega_{i,j,k} \int_{\Delta \Omega_{i,j}} d\ell d b \rho(s,\ell,b) s^2 \int_{4 \pi s^2 S_k^\text{min}}^{4 \pi s^2 S_k^\text{max}} {dN \over dL} dL \,,
$$ 

where $\omega_{i,j,k}$ is the efficiency factor for pular-like PS detection in that bin.  Note that we have assumed that $\omega_{i,j,k}$ is constant within the bin, so that we may bring it outside of the integral.
 

In the code, for both the disk and the Bulge, we normalize both $\rho$ and $dN/dL$ such that the full integral over all of space and flux is equal to the total number of sources $N$ for that population.  The number of sources is one of the model parameters.  To work with this convention, we must know how to properly normalize $\rho$ and $dN/dL$ .

The normalization of $dN/dL$ is common to both the disk and Bulge.  Let 

$$
{dN \over dL} = {N_0 \over L^\beta}
$$

for $L$ in the range $[L_\text{min}, L_\text{max}]$ and be zero outside of this range.  We will normalize $dN/dL$ so that it integrates to unity, which implies that
 
$$
N_0 = {\beta - 1 \over L_\text{min}^{1-\beta} - L_\text{max}^{1 - \beta} } \,.
$$  


We discuss $\rho$ now for the disk and Bulge.  We will use the normalization that $\int d^3x \rho = N$.

### Bulge Sources 

We characterize the Bulge PS population by

$$
\rho_\text{bulge} = {C \over r^\alpha} \,,
$$

where $r$ is the distance from the Galactic Center, so long as $r < r_\text{cut}$ .  $\rho$ is equal to zero for larger values of $r$ .  To have to proper normalization, it is straightforward to verify that
 
$$
C = {(3 - \alpha) N_\text{bulge} \over 4 \pi r_\text{cut}^{3 - \alpha} } \,.
$$


### Disk Sources

The disk is modeled by the distribution 

$$
\rho_\text{disk} = C R^n e^{-R / \sigma} e^{-|z| / z_0} \,,
$$

and we find that in this case

$$
C = {N_\text{disk} \over 4 \pi z_0 \sigma^{n+2} \Gamma(n+2)} \,.
$$


### Code Implementation

Here, we describe in more detail how we implement the computations of the expected number of PS counts in the file `get_counts_inline.pyx`.

The cartesian grid is given by 12 linear-spaced bins in the inner 40$^\circ$ of the Milky Way:
```python
cdef double[::1] ang_boundaries = np.linspace(-20.,20.,13,dtype=DTYPE)
``` 





 
