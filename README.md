# GCE-2FIG
This code uses 2FIG sources to constrain spatial distributions of pulsar-like PS populations.  The main purpose of the code is to test the results presented in 1705.00009.

## Compiling and Running

This code package is written in `python` and `cython`.  To compile the `cython`, go into the `python/` subfolder and execute the file `make.sh`.  Note that you must have `cython` installed.  Examples of how to run the code are presented in the `examples/` subfolder.

## The likelihood Function

The `GCE-2FIG` code package uses the likelihood function presented in 1705.00009.  The likelihood is imlemented in the file `likelihood.pyx` in the `python/` subfolder.  The sky is spatially binned in a cartesian grid, and the data consists of the number of pulsar candidates detected in each bin.  The model parameters characterize the disk and Bulge point source populations, and as we scan over the model parameters we calculated the expected number of detected sources in each bin, utlizing the _Fermi_-provided efficiency function for detecting sources at a given flux and spatial position.  The likelihood is then given by the product over all pixels of the Poisson probabilities to observe the data counts given the predicted model counts.

As an option, we also incorporate the prior distribution in 1705.00009.

## The Expected Number PS Counts

Here, we describe how we calculate the expected number of PS counts, given the model parameters.  These computations are implemented in the file `get_counts_inline.pyx` in the `python/` subfolder.  First, we describe the computation for a general PS population, with density function <img alt="$\rho(s,\ell,b)$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e33a91659757dc5a5a317888c76bc940.svg?c564d04929&invert_in_darkmode" align=middle width="57.31143pt" height="24.56553pt"/>, where <img alt="$(s,\ell,b)$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e4c7375533c0164f8ecc75cfe4198ea0.svg?a08062596&invert_in_darkmode" align=middle width="48.84429pt" height="24.56553pt"/> are spatial coordinates with origin at the Earth, and luminosity function <img alt="$dN/dL$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?bd96cccd52&invert_in_darkmode" align=middle width="49.96266pt" height="24.56553pt"/>.  Note that <img alt="$\ell$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d30a65b936d8007addc9c789d5a7ae49.svg?58e2d9969f&invert_in_darkmode" align=middle width="6.8238225pt" height="22.74591pt"/> and <img alt="$b$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?97c3419f40&invert_in_darkmode" align=middle width="7.0284885pt" height="22.74591pt"/> are Galactic longitude and latitude, respectively, while <img alt="$s$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?ce14f5a726&invert_in_darkmode" align=middle width="7.6767405pt" height="14.10255pt"/> is the distance from the Earth.  The expected number of counts in a given spatial bin (labeled by <img alt="$i,j$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4fe48dde86ac2d37419f0b35d57ac460.svg?8899b4f851&invert_in_darkmode" align=middle width="20.612625pt" height="21.60213pt"/>) and flux bin (labelled by <img alt="$k$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/63bb9849783d01d91403bc9a5fea12a2.svg?264269a38a&invert_in_darkmode" align=middle width="9.041505pt" height="22.74591pt"/>) is given by the integral

<p align="center"><img alt="$$&#10;N_{i,j,k} = \omega_{i,j,k} \int_{\Delta \Omega_{i,j}} d\ell d b \rho(s,\ell,b) s^2 \int_{4 \pi s^2 S_k^\text{min}}^{4 \pi s^2 S_k^\text{max}} {dN \over dL} dL \,,&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/0b66f12996fe020e588077671462ab65.svg?c6955b30d&invert_in_darkmode" align=middle width="384.76845pt" height="47.505645pt"/></p> 

where <img alt="$\omega_{i,j,k}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/2e71a536af409f25b2c3baebdba40859.svg?1f0dae30bf&invert_in_darkmode" align=middle width="35.22189pt" height="14.10255pt"/> is the efficiency factor for pular-like PS detection in that bin.  Note that we have assumed that <img alt="$\omega_{i,j,k}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/2e71a536af409f25b2c3baebdba40859.svg?d14a467875&invert_in_darkmode" align=middle width="35.22189pt" height="14.10255pt"/> is constant within the bin, so that we may bring it outside of the integral.
 

In the code, for both the disk and the Bulge, we normalize both <img alt="$\rho$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?8b6fc0573d&invert_in_darkmode" align=middle width="8.46714pt" height="14.10255pt"/> and <img alt="$dN/dL$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?680dd27ef4&invert_in_darkmode" align=middle width="49.96266pt" height="24.56553pt"/> such that the full integral over all of space and flux is equal to the total number of sources <img alt="$N$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f9c4988898e7f532b9f826a75014ed3c.svg?68afb844a6&invert_in_darkmode" align=middle width="14.94405pt" height="22.38192pt"/> for that population.  The number of sources is one of the model parameters.  To work with this convention, we must know how to properly normalize <img alt="$\rho$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?4c5784669b&invert_in_darkmode" align=middle width="8.46714pt" height="14.10255pt"/> and <img alt="$dN/dL$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?87aadfc47d&invert_in_darkmode" align=middle width="49.96266pt" height="24.56553pt"/> .

The normalization of <img alt="$dN/dL$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?4f13d55fe5&invert_in_darkmode" align=middle width="49.96266pt" height="24.56553pt"/> is common to both the disk and Bulge.  Let 

<p align="center"><img alt="$$&#10;{dN \over dL} = {N_0 \over L^\beta}&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/dc7e9cac24bdd8b9fb9aa0389e209134.svg?38974d1a7b&invert_in_darkmode" align=middle width="69.953235pt" height="33.769395pt"/></p>

for <img alt="$L$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ddcb483302ed36a59286424aa5e0be17.svg?959f71e49&invert_in_darkmode" align=middle width="11.14542pt" height="22.38192pt"/> in the range <img alt="$[L_\text{min}, L_\text{max}]$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f8702ab460607cfc76e52498178edb88.svg?51c782a70e&invert_in_darkmode" align=middle width="86.361pt" height="24.56553pt"/> and be zero outside of this range.  We will normalize <img alt="$dN/dL$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?256ac45d0&invert_in_darkmode" align=middle width="49.96266pt" height="24.56553pt"/> so that it integrates to unity, which implies that
 
<p align="center"><img alt="$$&#10;N_0 = {\beta - 1 \over L_\text{min}^{1-\beta} - L_\text{max}^{1 - \beta} } \,.&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/0a8a3c1792e9fdb071d6a6920dec93d7.svg?1ef054e12e&invert_in_darkmode" align=middle width="147.489705pt" height="41.283165pt"/></p>  


We discuss <img alt="$\rho$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?c766a6f590&invert_in_darkmode" align=middle width="8.46714pt" height="14.10255pt"/> now for the disk and Bulge.  We will use the normalization that <img alt="$\int d^3x \rho = N$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/91f9537316cdc8858af8eeb325936c77.svg?88afd2edfa&invert_in_darkmode" align=middle width="84.237945pt" height="26.70657pt"/>.

### Bulge Sources 

We characterize the Bulge PS population by

<p align="center"><img alt="$$&#10;\rho_\text{bulge} = {C \over r^\alpha} \,,&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/20a7cb8590592ee1ab7830a4d8541614.svg?3fe74a9190&invert_in_darkmode" align=middle width="90.27315pt" height="33.5874pt"/></p>

where <img alt="$r$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/89f2e0d2d24bcf44db73aab8fc03252c.svg?78a8d41a23&invert_in_darkmode" align=middle width="7.8435885pt" height="14.10255pt"/> is the distance from the Galactic Center, so long as <img alt="$r &lt; r_\text{cut}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/af61c60decc49acae8eeda9dca8fe898.svg?2fbe55376f&invert_in_darkmode" align=middle width="55.277805pt" height="17.65764pt"/> .  <img alt="$\rho$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?670db60887&invert_in_darkmode" align=middle width="8.46714pt" height="14.10255pt"/> is equal to zero for larger values of <img alt="$r$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/89f2e0d2d24bcf44db73aab8fc03252c.svg?ab9525e33c&invert_in_darkmode" align=middle width="7.8435885pt" height="14.10255pt"/> .  To have to proper normalization, it is straightforward to verify that
 
<p align="center"><img alt="$$&#10;C = {(3 - \alpha) N_\text{bulge} \over 4 \pi r_\text{cut}^{3 - \alpha} } \,.&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/9e79c22983590f5b42a2672ec8cf10a1.svg?aa4827ab3d&invert_in_darkmode" align=middle width="142.32603pt" height="40.08807pt"/></p>


### Disk Sources

The disk is modeled by the distribution 

<p align="center"><img alt="$$&#10;\rho_\text{disk} = C R^n e^{-R / \sigma} e^{-|z| / z_0} \,,&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ffd30829e06c1a1f6e2acec3a70ed61c.svg?c5b702842f&invert_in_darkmode" align=middle width="192.0963pt" height="18.569595pt"/></p>

and we find that in this case

<p align="center"><img alt="$$&#10;C = {N_\text{disk} \over 4 \pi z_0 \sigma^{n+2} \Gamma(n+2)} \,.&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/a5bba1689089f869eba4d9fda47f21e7.svg?dc407a75c5&invert_in_darkmode" align=middle width="175.99395pt" height="37.68171pt"/></p>


### Code Implementation

Here, we describe in more detail how we implement the computations of the expected number of PS counts in the file `get_counts_inline.pyx`.

The cartesian grid is given by 12 linear-spaced bins in the inner 40<img alt="$^\circ$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bda93e7eec1ea3bd03d7177c5b991481.svg?c706b6e261&invert_in_darkmode" align=middle width="6.7100715pt" height="22.59873pt"/> of the Milky Way:

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

return the number of counts in the Bulge in spatial <img alt="$\ell$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d30a65b936d8007addc9c789d5a7ae49.svg?908fa2bcc9&invert_in_darkmode" align=middle width="6.8238225pt" height="22.74591pt"/> bin `i`, <img alt="$b$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?8e4c5b8dbc&invert_in_darkmode" align=middle width="7.0284885pt" height="22.74591pt"/> bin `j`, and flux bin `k`.  Most of the parameters above are defined previously.  However, there are a few parameters that are specific to the calculation method.  These include

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

Then, the bulk of the code proceeds through the loops used to perform the angular integration and the integration over <img alt="$s$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?b776b1e10a&invert_in_darkmode" align=middle width="7.6767405pt" height="14.10255pt"/>:

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

Above, `i_ell` is the index over longitude, `i_b` is the index of latitude, and `l` is the index over the <img alt="$s$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?61cd46dfa9&invert_in_darkmode" align=middle width="7.6767405pt" height="14.10255pt"/> bins.  In the first few lines we simply update the new values for `b` and `ell` and their cosines.  We then see if these values of `b` and `ell` are masked.  The statment `if incl > 0:` looks to see if we are too far away from the GC, such that we will completely miss the bulge.  This is relevant given the finite extent of the bulge, but this step is not needed for the disk contribution.  Then, we begin looping over the `s` integral, through the index `l`.  First, we find the true limits of the flux integral, given that fact that the luminosity function is cut off at `Lmin` and `Lmax`.  Then, we perform the integral over `L` to get `pref_L`.  The quantity `integral` contains both the integral over `s` and the integral over `L`.  

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







 
