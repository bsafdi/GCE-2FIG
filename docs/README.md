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

```c
cdef double[::1] ang_boundaries = np.linspace(-20.,20.,13,dtype=DTYPE)
``` 

The flux-bin boundaries are given by

```python
fluxvals = [1.00000000e-06, 1.46779927e-06, 2.15443469e-06, 3.16227766e-06, 4.64158883e-06, 6.81292069e-06, 1.00000000e-05, 3.16227766e-05, 1.00000000e-04]
```

in units of MeV/cm^2/s.

The function 

```c
double Nbulge_full_ang_ijk(int i, int j, int k, double Nbulge, double omega_ijk, double alpha, double beta,double rcut , double Lmin, double Lmax ,int Ns ,int Nang, double theta_mask )
```

return the number of counts in the Bulge in spatial $\ell$ bin `i`, $b$ bin `j`, and flux bin `k`.  Most of the parameters above are defined previously.  However, there are a few parameters that are specific to the calculation method.  These include

    1. int Ns: the number of points in the numerical integral over s
    2. int Nang: the number of latitude and longitude points in the numerical integral over the pixel
    3. double theta_mask: the radius in degrees fro the Galactic Center that will be masked when computing the number of counts
  
The code itself proceeds by first initializing relevant quantities, such as

```python
cdef double ell_start = ang_boundaries[i]
cdef double b_start = ang_boundaries[j]
cdef double d_ang = (ang_boundaries[1] - ang_boundaries[0])/float(Nang)
```

and  calculating relevant pre-factors, such as 

```python
pref_rho = omega_ijk * Nbulge * (3. - alpha) / 4. / pi / pow(rcut,3 - alpha)
```

and 

```python
pref_L_norm = 1./(pow(Lmin,1-beta) - pow(Lmax,1-beta))
```

Then, the bulk of the code proceeds through the loops used to perform the angular integration and the integration over $s$:

```python
for i_b in range(Nang):
    b = b_start + i_b*d_ang + d_ang/2.
    for i_ell in range(Nang):
        ell = ell_start + i_ell*d_ang + d_ang/2.
        coslval = cos(ell * degtorad)
        cosbval = cos(b * degtorad)

        if cosbval * coslval < costhetaval: #is it masked?

            incl = 1. - (1. - pow(rcut/rodot,2) ) * pow(coslval*cosbval,2)
            if incl > 0: # should we bother, or is the answer going to be zero?
                integral = 0.0
                s = smin
                for l in range(Ns): #Loop over `s` for `s` integral

                    flux_min_eval = max(Lmin/s**2,4*pi*fluxmin)
                    flux_max_eval = min(Lmax/s**2,4*pi*fluxmax)

                    if flux_max_eval < flux_min_eval:
                        pref_L = 0.0
                    else:
                        pref_L = (pow(flux_min_eval, 1.-beta) - pow(flux_max_eval, 1.-beta))*pref_L_norm

                    r_squared = (s*cosbval*coslval - rodot)**2 + s**2*(1. - pow(cosbval*coslval,2) )
                    if r_squared < pow(rcut,2):
                        integral += pref_L * pow(s,4.-2.*beta)*pow(r_squared,-alpha/2.)
                    s += ds

                total_res += cosbval * integral 
```

Above, `i_ell` is the index over longitude, `i_b` is the index of latitude, and `l` is the index over the $s$ bins.  In the first few lines we simply update the new values for `b` and `ell` and their cosines.  We then see if these values of `b` and `ell` are masked.  The statment `if incl > 0:` looks to see if we are too far away from the GC, such that we will completely miss the bulge.  This is relevant given the finite extent of the bulge, but this step is not needed for the disk contribution.  Then, we begin looping over the `s` integral, through the index `l`.  First, we find the true limits of the flux integral, given that fact that the luminosity function is cut off at `Lmin` and `Lmax`.  Then, we perform the integral over `L` to get `pref_L`.  The quantity `integral` contains both the integral over `s` and the integral over `L`.  

Finally, we return 

```python
return total_res * d_ang**2. * degtorad**2 * ds * pref_rho
```  

Note that the disk proceeds similarly, as do the functions that compute the full-sky integrals for the prior distribution.
    
## Compiling this document

To render this on `github`, execute the following:
```
python -m readme2tex --rerender --bustcache --output README.md docs/README.md
```

# Results

Unless floated, $n = 2.35$,
$\sigma = 1.528$,
$\alpha = 2.6$,
$\beta = 1.2$,
$z_0 = 0.7$.

Only floating $N_D$ and $N_B$:

| $N_D$        | $N_B$       | TS   |
|--------------|-------------|------|
| $4462\pm635$ | 0           | 0    |
| $3152\pm800$ | $690\pm359$ | 9    |

Floating other paramsters as in 2FIG Table 2:


| $N_D$             | $N_B$                     | $z_0$         | $\beta$       | $\alpha$      | TS |
|-------------------|---------------------------|---------------|---------------|---------------|----|
| $953805\pm149704$ | 0                         | $0.11\pm0.05$ | $2\pm0.08$    | -             | 0  |
| $675207\pm174537$ | $265164\pm117051$         | $0.07\pm0.04$ | $2\pm0.14$    | -             | 11 |
| $869369\pm426897$ | $849551\pm2.53\times10^6$ | $0.08\pm0.06$ | $1.99\pm0.72$ | $2.93\pm0.66$ | 15 |








 
