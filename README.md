# GCE-2FIG
This code uses 2FIG sources to constrain spatial distributions of pulsar-like PS populations.  The main purpose of the code is to test the results presented in [1705.00009](https://arxiv.org/pdf/1705.00009.pdf).

## Compiling and running

This code package is written in `python` and `cython`.  To compile the `cython`, go into the `likelihood/` subfolder and execute the file `make.sh`.  Note that you must have `cython` installed.  Examples of how to run the code are presented in the `examples/` subfolder.

## The data

## The likelihood function

The `GCE-2FIG` code package uses the likelihood function presented in [1705.00009](https://arxiv.org/pdf/1705.00009.pdf).  The likelihood is implemented in the file `likelihood.pyx` in the `likelihood/` subfolder.  The sky is spatially binned in a cartesian grid, and the data consists of the number of pulsar candidates detected in each bin.  The model parameters characterize the disk and Bulge point source populations, and as we scan over the model parameters we calculated the expected number of detected sources in each bin, utilizing the _Fermi_-provided efficiency function for detecting sources at a given flux and spatial position.  The likelihood is then given by the product over all pixels of the Poisson probabilities to observe the data counts given the predicted model counts.

As an option, we also incorporate the prior distribution in [1705.00009](https://arxiv.org/pdf/1705.00009.pdf).

## The expected number of PS counts

Here, we describe how we calculate the expected number of PS counts, given the model parameters.  These computations are implemented in the file `get_counts_inline.pyx` in the `likelihood/` subfolder.  First, we describe the computation for a general PS population, with density function <img alt="$\rho(s,\ell,b)$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e33a91659757dc5a5a317888c76bc940.svg?135662e0b0&invert_in_darkmode" align=middle width="57.31143pt" height="24.56553pt"/>, where <img alt="$(s,\ell,b)$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e4c7375533c0164f8ecc75cfe4198ea0.svg?7fa93e0114&invert_in_darkmode" align=middle width="48.84429pt" height="24.56553pt"/> are spatial coordinates with origin at the Earth, and luminosity function <img alt="$dN/dL$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?394f58780d&invert_in_darkmode" align=middle width="49.96266pt" height="24.56553pt"/>.  Note that <img alt="$\ell$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d30a65b936d8007addc9c789d5a7ae49.svg?1cf828d251&invert_in_darkmode" align=middle width="6.8238225pt" height="22.74591pt"/> and <img alt="$b$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?3995720df4&invert_in_darkmode" align=middle width="7.0284885pt" height="22.74591pt"/> are Galactic longitude and latitude, respectively, while <img alt="$s$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?cf818cc893&invert_in_darkmode" align=middle width="7.6767405pt" height="14.10255pt"/> is the distance from the Earth.  The expected number of counts in a given spatial bin (labeled by <img alt="$i,j$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4fe48dde86ac2d37419f0b35d57ac460.svg?c87608b475&invert_in_darkmode" align=middle width="20.612625pt" height="21.60213pt"/>) and flux bin (labelled by <img alt="$k$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/63bb9849783d01d91403bc9a5fea12a2.svg?6ca70f7740&invert_in_darkmode" align=middle width="9.041505pt" height="22.74591pt"/>) is given by the integral

<p align="center"><img alt="$$&#10;N_{i,j,k} = \omega_{i,j,k} \int_{\Delta \Omega_{i,j}} d\ell d b \rho(s,\ell,b) s^2 \int_{4 \pi s^2 S_k^\text{min}}^{4 \pi s^2 S_k^\text{max}} {dN \over dL} dL \,,&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/0b66f12996fe020e588077671462ab65.svg?9b6ad4bf19&invert_in_darkmode" align=middle width="384.76845pt" height="47.505645pt"/></p> 

where <img alt="$\omega_{i,j,k}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/2e71a536af409f25b2c3baebdba40859.svg?8076a65aac&invert_in_darkmode" align=middle width="35.22189pt" height="14.10255pt"/> is the efficiency factor for pulsar-like PS detection in that bin.  Note that we have assumed that <img alt="$\omega_{i,j,k}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/2e71a536af409f25b2c3baebdba40859.svg?d353eadf58&invert_in_darkmode" align=middle width="35.22189pt" height="14.10255pt"/> is constant within the bin, so that we may bring it outside of the integral.
 

In the code, for both the disk and the Bulge, we normalize both <img alt="$\rho$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?aa7e3b606e&invert_in_darkmode" align=middle width="8.46714pt" height="14.10255pt"/> and <img alt="$dN/dL$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?551676c132&invert_in_darkmode" align=middle width="49.96266pt" height="24.56553pt"/> such that the full integral over all of space and flux is equal to the total number of sources <img alt="$N$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f9c4988898e7f532b9f826a75014ed3c.svg?d76d5be907&invert_in_darkmode" align=middle width="14.94405pt" height="22.38192pt"/> for that population.  The number of sources is one of the model parameters.  To work with this convention, we must know how to properly normalize <img alt="$\rho$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?81a0f149b8&invert_in_darkmode" align=middle width="8.46714pt" height="14.10255pt"/> and <img alt="$dN/dL$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?5eb20a56cc&invert_in_darkmode" align=middle width="49.96266pt" height="24.56553pt"/> .

The normalization of <img alt="$dN/dL$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?19cb08d2c7&invert_in_darkmode" align=middle width="49.96266pt" height="24.56553pt"/> is common to both the disk and Bulge.  Let 

<p align="center"><img alt="$$&#10;{dN \over dL} = {N_0 \over L^\beta}&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/dc7e9cac24bdd8b9fb9aa0389e209134.svg?91938cbc29&invert_in_darkmode" align=middle width="69.953235pt" height="33.769395pt"/></p>

for <img alt="$L$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ddcb483302ed36a59286424aa5e0be17.svg?64a2874f6a&invert_in_darkmode" align=middle width="11.14542pt" height="22.38192pt"/> in the range <img alt="$[L_\text{min}, L_\text{max}]$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f8702ab460607cfc76e52498178edb88.svg?46c9d8ecac&invert_in_darkmode" align=middle width="86.361pt" height="24.56553pt"/> and be zero outside of this range.  We will normalize <img alt="$dN/dL$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?304be09304&invert_in_darkmode" align=middle width="49.96266pt" height="24.56553pt"/> so that it integrates to unity, which implies that
 
<p align="center"><img alt="$$&#10;N_0 = {\beta - 1 \over L_\text{min}^{1-\beta} - L_\text{max}^{1 - \beta} } \,.&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/0a8a3c1792e9fdb071d6a6920dec93d7.svg?6a971a0633&invert_in_darkmode" align=middle width="147.489705pt" height="41.283165pt"/></p>  


We discuss <img alt="$\rho$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?539e299e8&invert_in_darkmode" align=middle width="8.46714pt" height="14.10255pt"/> now for the disk and Bulge.  We will use the normalization that <img alt="$\int d^3x \rho = N$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/91f9537316cdc8858af8eeb325936c77.svg?14ea8a958e&invert_in_darkmode" align=middle width="84.237945pt" height="26.70657pt"/>.

### Bulge sources 

We characterize the Bulge PS population by

<p align="center"><img alt="$$&#10;\rho_\text{bulge} = {C \over r^\alpha} \,,&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/20a7cb8590592ee1ab7830a4d8541614.svg?41716d1293&invert_in_darkmode" align=middle width="90.27315pt" height="33.5874pt"/></p>

where <img alt="$r$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/89f2e0d2d24bcf44db73aab8fc03252c.svg?90f3d05a14&invert_in_darkmode" align=middle width="7.8435885pt" height="14.10255pt"/> is the distance from the Galactic Center, so long as <img alt="$r &lt; r_\text{cut}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/af61c60decc49acae8eeda9dca8fe898.svg?11b3ab743b&invert_in_darkmode" align=middle width="55.277805pt" height="17.65764pt"/> .  <img alt="$\rho$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?bb6e9c02b8&invert_in_darkmode" align=middle width="8.46714pt" height="14.10255pt"/> is equal to zero for larger values of <img alt="$r$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/89f2e0d2d24bcf44db73aab8fc03252c.svg?ad44d92f6c&invert_in_darkmode" align=middle width="7.8435885pt" height="14.10255pt"/> .  To have to proper normalization, it is straightforward to verify that
 
<p align="center"><img alt="$$&#10;C = {(3 - \alpha) N_\text{bulge} \over 4 \pi r_\text{cut}^{3 - \alpha} } \,.&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/9e79c22983590f5b42a2672ec8cf10a1.svg?e4f528e4cd&invert_in_darkmode" align=middle width="142.32603pt" height="40.08807pt"/></p>


### Disk sources

The disk is modeled by the distribution 

<p align="center"><img alt="$$&#10;\rho_\text{disk} = C R^n e^{-R / \sigma} e^{-|z| / z_0} \,,&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ffd30829e06c1a1f6e2acec3a70ed61c.svg?b0025cc5be&invert_in_darkmode" align=middle width="192.0963pt" height="18.569595pt"/></p>

and we find that in this case

<p align="center"><img alt="$$&#10;C = {N_\text{disk} \over 4 \pi z_0 \sigma^{n+2} \Gamma(n+2)} \,.&#10;$$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/a5bba1689089f869eba4d9fda47f21e7.svg?a14ae9f465&invert_in_darkmode" align=middle width="175.99395pt" height="37.68171pt"/></p>


### Code implementation

Here, we describe in more detail how we implement the computations of the expected number of PS counts in the file `get_counts_inline.pyx`.

The cartesian grid is given by 12 linear-spaced bins in the inner 40<img alt="$^\circ$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bda93e7eec1ea3bd03d7177c5b991481.svg?8e3736ea45&invert_in_darkmode" align=middle width="6.7100715pt" height="22.59873pt"/> of the Milky Way:

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

return the number of counts in the Bulge in spatial <img alt="$\ell$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d30a65b936d8007addc9c789d5a7ae49.svg?d6c0dce3a8&invert_in_darkmode" align=middle width="6.8238225pt" height="22.74591pt"/> bin `i`, <img alt="$b$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?86633de76b&invert_in_darkmode" align=middle width="7.0284885pt" height="22.74591pt"/> bin `j`, and flux bin `k`.  Most of the parameters above are defined previously.  However, there are a few parameters that are specific to the calculation method.  These include

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

Then, the bulk of the code proceeds through the loops used to perform the angular integration and the integration over <img alt="$s$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?9ea41d6366&invert_in_darkmode" align=middle width="7.6767405pt" height="14.10255pt"/>:

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

Above, `i_ell` is the index over longitude, `i_b` is the index of latitude, and `l` is the index over the <img alt="$s$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?51c50973df&invert_in_darkmode" align=middle width="7.6767405pt" height="14.10255pt"/> bins.  In the first few lines we simply update the new values for `b` and `ell` and their cosines.  We then see if these values of `b` and `ell` are masked.  The statement `if incl > 0:` looks to see if we are too far away from the GC, such that we will completely miss the bulge.  This is relevant given the finite extent of the bulge, but this step is not needed for the disk contribution.  Then, we begin looping over the `s` integral, through the index `l`.  First, we find the true limits of the flux integral, given that fact that the luminosity function is cut off at `Lmin` and `Lmax`.  Then, we perform the integral over `L` to get `pref_L`.  The quantity `integral` contains both the integral over `s` and the integral over `L`.  

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

## Results

Here, we present some of the main results from the analyses performed with this code package.  We will use three different efficiency functions.

	1.  A longitude independent efficiency function, supplied directly in [[1705.00009](https://arxiv.org/pdf/1705.00009.pdf)](https://arxiv.org/pdf/[1705.00009](https://arxiv.org/pdf/1705.00009.pdf).pdf).  We will label this as `long-indep`
	2. A longitude-corrected efficiency function, corrected using data in [1305.4385](https://arxiv.org/abs/1305.4385).  We will label this as (`2PC-corrected`)
	3. A longitude-corrected efficiency function, corrected by rescaling the longitude-independent efficiency function by the actual data values within the pixels.  We will label this as (`data-corrected`)

The reason we use three different efficiency functions is that this allows us to understand the systematic uncertainty associated with the fact that we likely do not have the correct efficiency.

Below, unless floated, values are taken to be:

- <img alt="$n = 2.35$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8d355ff26fea5f8b634222b72f19472f.svg?cfdb825bce&invert_in_darkmode" align=middle width="60.814545pt" height="21.10812pt"/>,
- <img alt="$\sigma = 1.528$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/5a341cc43b8429c286e732c54f035ae3.svg?9989a33b59&invert_in_darkmode" align=middle width="69.12081pt" height="21.10812pt"/>,
- <img alt="$\alpha = 2.6$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/fa32ade5e29240b938ede1a9fe79cd17.svg?ab8d8c130&invert_in_darkmode" align=middle width="53.33328pt" height="21.10812pt"/>,
- <img alt="$\beta = 1.2$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/01a4d380c848b0b24c9d376c25ec7880.svg?26901165b3&invert_in_darkmode" align=middle width="52.926885pt" height="22.74591pt"/>,
- <img alt="$z_0 = 0.7$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f2ca6be66d1c9df73aac3a892fda5efe.svg?17068737d7&invert_in_darkmode" align=middle width="57.78663pt" height="21.10812pt"/>.

Priors used: 

- <img alt="$N_B, N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/baca37ec0463c58e9b4bca0c53e033d2.svg?71e3fdbae2&invert_in_darkmode" align=middle width="55.98087pt" height="22.38192pt"/>: [0,3000000]
- <img alt="$\alpha$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c745b9b57c145ec5577b82542b2df546.svg?7aac4fd9ff&invert_in_darkmode" align=middle width="10.537065pt" height="14.10255pt"/>: [2.1, 5.0]
- <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?a666b0ad35&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/>: [0.01, 2.0]
- <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?60788f9786&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/>: [1.1,3.0]

Floating parameters as in Table 2 of [[1705.00009](https://arxiv.org/pdf/1705.00009.pdf)](https://arxiv.org/pdf/[1705.00009](https://arxiv.org/pdf/1705.00009.pdf).pdf):

<!--### Minuit

| <img alt="$N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?42e57d2870&invert_in_darkmode" align=middle width="24.21903pt" height="22.38192pt"/>             | <img alt="$N_B$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?41f74d2f15&invert_in_darkmode" align=middle width="23.61183pt" height="22.38192pt"/>                     | <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?6da2775a9b&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/>         | <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?8dcde15503&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/>       | <img alt="$\alpha$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c745b9b57c145ec5577b82542b2df546.svg?26733f8fa0&invert_in_darkmode" align=middle width="10.537065pt" height="14.10255pt"/>      | TS |
|-------------------|---------------------------|---------------|---------------|---------------|----|
| <img alt="$(1.53\pm0.76)\times 10^6$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/56708071c8a2440d9772b0df0b7504b6.svg?80a889a4c7&invert_in_darkmode" align=middle width="133.959045pt" height="26.70657pt"/> | <img alt="$0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/29632a9bf827ce0200454dd32fc3be82.svg?438e61ea0c&invert_in_darkmode" align=middle width="8.188554pt" height="21.10812pt"/>                         | <img alt="$0.06\pm0.05$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e8cc27253cf12d74d0ea54c7bd3e158a.svg?70d65f1023&invert_in_darkmode" align=middle width="78.272865pt" height="21.10812pt"/> | <img alt="$2.10\pm0.11$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/a4bb9a4c769c2c81d1565bf72c22d4ae.svg?5a83d26e71&invert_in_darkmode" align=middle width="78.272865pt" height="21.10812pt"/>    | -             | <img alt="$0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/29632a9bf827ce0200454dd32fc3be82.svg?e4ef1241ab&invert_in_darkmode" align=middle width="8.188554pt" height="21.10812pt"/>  |
| <img alt="$(9.00\pm5.24)\times 10^5$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/debb7cc9ede370461e0470142cd3fdf4.svg?c81c22bb84&invert_in_darkmode" align=middle width="133.959045pt" height="26.70657pt"/>    | <img alt="$(3.76\pm3.24)\times 10^5$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/a0403e487e8f755bbe24de662710518b.svg?cbe346023&invert_in_darkmode" align=middle width="133.959045pt" height="26.70657pt"/>      | <img alt="$ 0.05\pm0.04$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c894042905bb67272da0d1d501cad8aa.svg?6a31b39201&invert_in_darkmode" align=middle width="78.272865pt" height="21.10812pt"/>| <img alt="$ 2.05\pm0.10$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/91ddfe0948bd73958c7548e2c59450dc.svg?5df1861299&invert_in_darkmode" align=middle width="78.272865pt" height="21.10812pt"/>    | <img alt="$2.6$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/7509c44aea1098662d4f444907d5a427.svg?4815430b86&invert_in_darkmode" align=middle width="20.92629pt" height="21.10812pt"/>             | <img alt="$10.37$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/639548d548b582f9c90c4577fde9ef19.svg?e3a7a97c15&invert_in_darkmode" align=middle width="37.3032pt" height="21.10812pt"/> |
| <img alt="$(9.17\pm4.98)\times 10^5$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/5a2b0879440cf60d11da89ea650918f7.svg?7235795d81&invert_in_darkmode" align=middle width="133.959045pt" height="26.70657pt"/> | <img alt="$(2.99\pm1.53)\times 10^5$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/51e9145419309c59bb0cd8e8efc89d1e.svg?ba0899d839&invert_in_darkmode" align=middle width="133.959045pt" height="26.70657pt"/> | <img alt="$0.05\pm0.04$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8e1d270a78c11ca743bbcb8466b806e5.svg?a6aae34360&invert_in_darkmode" align=middle width="78.272865pt" height="21.10812pt"/> | <img alt="$2.06\pm0.10$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/9a035f597e130bd1552ca529d4ae5128.svg?184f5b9dd1&invert_in_darkmode" align=middle width="78.272865pt" height="21.10812pt"/> | <img alt="$2.96\pm0.03$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/02f516266e43f9ddc0827472462b4641.svg?48fa10524d&invert_in_darkmode" align=middle width="78.272865pt" height="21.10812pt"/> | <img alt="$11.66$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e80087059c9607f2cd892a1a1fb74283.svg?5724d97d68&invert_in_darkmode" align=middle width="37.3032pt" height="21.10812pt"/> |
-->

All results below are computed using `multinest`.  
### `Long-Indep`

| <img alt="$N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?1ca8b84c88&invert_in_darkmode" align=middle width="24.21903pt" height="22.38192pt"/>             | <img alt="$N_B$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?24b8d32764&invert_in_darkmode" align=middle width="23.61183pt" height="22.38192pt"/>                     | <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?b37309667f&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/>         | <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?c2c529e4ab&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/>       | <img alt="$\alpha$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c745b9b57c145ec5577b82542b2df546.svg?27f59be431&invert_in_darkmode" align=middle width="10.537065pt" height="14.10255pt"/>      | TS |
|-------------------|---------------------------|---------------|---------------|---------------|----|
| <img alt="$1483854^{+555529}_{-453093}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c2abccecbf075879688fff443609a78a.svg?67760e7d03&invert_in_darkmode" align=middle width="106.72365pt" height="28.83969pt"/> | <img alt="$0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/29632a9bf827ce0200454dd32fc3be82.svg?39ca0a1429&invert_in_darkmode" align=middle width="8.188554pt" height="21.10812pt"/>                         | <img alt="$0.08^{+0.05}_{-0.03}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/a6709f75746c9fe2431f875bc873fac5.svg?386e17617c&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/> | <img alt="$2.09^{+0.07}_{-0.08}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/2fc4006f8e81856f645d477be8c15c89.svg?70dacbd118&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/>    | -             | <img alt="$0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/29632a9bf827ce0200454dd32fc3be82.svg?60e53e9ade&invert_in_darkmode" align=middle width="8.188554pt" height="21.10812pt"/>  |
|<img alt="$1153102^{+525377}_{-378358}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/89ac6a656a5fe30c65937b0c85cf6430.svg?66315e9d71&invert_in_darkmode" align=middle width="106.72365pt" height="28.83969pt"/>      | <img alt="$506568^{+496453}_{-238196}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/b07986458776d26645734f823884e753.svg?aa18637ab7&invert_in_darkmode" align=middle width="98.53503pt" height="28.83969pt"/>    |  <img alt="$0.05^{+0.03}_{-0.02}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/31bc73fa691ca22a146f69461671e19e.svg?53f85eadfb&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/>| <img alt="$2.10^{+0.08}_{-0.07}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/9f056d70bba38553c351ba80f246f086.svg?4f541b28c6&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/>    | <img alt="$2.6$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/7509c44aea1098662d4f444907d5a427.svg?183d2e82ba&invert_in_darkmode" align=middle width="20.92629pt" height="21.10812pt"/>             | <img alt="$10.34$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/aac197995b881999d6d8a8cb8ae60c5a.svg?4508112f5&invert_in_darkmode" align=middle width="37.3032pt" height="21.10812pt"/> |
| <img alt="$1052708^{+456557}_{-343513}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/5943d75e21b33d29b9c03e0ddbe64211.svg?b3fc1605b2&invert_in_darkmode" align=middle width="106.72365pt" height="28.83969pt"/> | <img alt="$966064^{+1160994}_{-587532}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/dfbfd8205ddc44e1982a78fc4c668aed.svg?312a192ebd&invert_in_darkmode" align=middle width="104.88126pt" height="28.83969pt"/> | <img alt="$0.06^{+0.04}_{-0.02}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/5eaf3c8b3988b60f7568e7a20f266d69.svg?90142e6b37&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/> | <img alt="$2.08^{+0.07}_{-0.07}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ce577a0d363bb29ad852707f18305b22.svg?26474ed4a7&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/> | <img alt="$2.83^{+0.11}_{-0.27}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/97014f8146a8320e7598bfabad3ddbf7.svg?97a36a52d8&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/> | <img alt="$11.39$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/983da5045b8036e65ea66cf00a6b72ae.svg?e0eb5c5aeb&invert_in_darkmode" align=middle width="37.3032pt" height="21.10812pt"/> |


#### Float <img alt="$N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?e622888390&invert_in_darkmode" align=middle width="24.21903pt" height="22.38192pt"/>, <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?129bbc29b2&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/> and <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?a787bba46e&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/>

![Disk only](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-z0-beta-lon-indep_eff.png "Disk only")

#### Float <img alt="$N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?598c30136d&invert_in_darkmode" align=middle width="24.21903pt" height="22.38192pt"/>, <img alt="$N_B$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?b833690ff&invert_in_darkmode" align=middle width="23.61183pt" height="22.38192pt"/>, <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?28b09a4752&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/> and <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?9539ee79ef&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/>

![Disk and bulge](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-nb-z0-beta-lon-indep_eff.png "Disk and bulge")

#### Float <img alt="$N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?8ad39ccf0&invert_in_darkmode" align=middle width="24.21903pt" height="22.38192pt"/>, <img alt="$N_B$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?2674caa35b&invert_in_darkmode" align=middle width="23.61183pt" height="22.38192pt"/>, <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?8e2dbb4819&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/>, <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?12c6c54ade&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/> and <img alt="$\alpha$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c745b9b57c145ec5577b82542b2df546.svg?e5e114784c&invert_in_darkmode" align=middle width="10.537065pt" height="14.10255pt"/>

![Float alpha](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-nb-alpha-z0-beta-lon-indep_eff.png "Float alpha")


### `2PC-Corrected`

| <img alt="$N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?a71ae2664d&invert_in_darkmode" align=middle width="24.21903pt" height="22.38192pt"/>             | <img alt="$N_B$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?e67a14e17&invert_in_darkmode" align=middle width="23.61183pt" height="22.38192pt"/>                     | <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?8e2f6ffc7e&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/>         | <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?973b354717&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/>       | <img alt="$\alpha$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c745b9b57c145ec5577b82542b2df546.svg?b4e2211c54&invert_in_darkmode" align=middle width="10.537065pt" height="14.10255pt"/>      | TS |
|-------------------|---------------------------|---------------|---------------|---------------|----|
| <img alt="$1482073^{+575204}_{-465966}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4430b22ef4ce0098b272dfb31846a981.svg?671a0e64e4&invert_in_darkmode" align=middle width="106.72365pt" height="28.83969pt"/> | <img alt="$0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/29632a9bf827ce0200454dd32fc3be82.svg?d5911aefd7&invert_in_darkmode" align=middle width="8.188554pt" height="21.10812pt"/>                         | <img alt="$0.08^{+0.05}_{-0.03}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/a6709f75746c9fe2431f875bc873fac5.svg?986615cfd8&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/> | <img alt="$2.09_{-0.08}^{+0.08}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/326a7cda247ae13fa0ee7fac8829a4dd.svg?c50f8aca72&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/>    | -             | 0  |
|<img alt="$1068731^{+483159}_{-355819}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/a3056fa9e7183875d5803231c4f25888.svg?591104e117&invert_in_darkmode" align=middle width="106.72365pt" height="28.83969pt"/>      | <img alt="$590663^{+540377}_{-282268}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6866b8b029f55318306668f70563aa0f.svg?59452eec86&invert_in_darkmode" align=middle width="98.53503pt" height="28.83969pt"/>    |  <img alt="$0.05^{+0.04}_{-0.02}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/7f2b9cc04db42ea293fbdf287e499878.svg?3b69fcc3c2&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/>| <img alt="$2.10^{+0.08}_{-0.08}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/38fe32d8bb8dee1079a1c351e829745e.svg?173b093f63&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/>    | <img alt="$2.6$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/7509c44aea1098662d4f444907d5a427.svg?6c4d7cf72f&invert_in_darkmode" align=middle width="20.92629pt" height="21.10812pt"/>             | <img alt="$13.21$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e7f87f6aaa9faee4ac076a9264bfa1b6.svg?723a5a774e&invert_in_darkmode" align=middle width="37.3032pt" height="21.10812pt"/> |
| <img alt="$955272^{+409682}_{-318262}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c5878c32063616c2d7a6678d4b3daa21.svg?896e0b7243&invert_in_darkmode" align=middle width="98.53503pt" height="28.83969pt"/> | <img alt="$1328045.71^{+1066159}_{-803946}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/45b7829e196c71fa160fc9499f1fe512.svg?716f40e195&invert_in_darkmode" align=middle width="133.99584pt" height="28.83969pt"/> | <img alt="$0.06^{+0.04}_{-0.02}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/5eaf3c8b3988b60f7568e7a20f266d69.svg?9af7d54c6d&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/> | <img alt="$2.08^{+0.07}_{-0.07}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ce577a0d363bb29ad852707f18305b22.svg?4250af670&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/> | <img alt="$2.86^{+0.08}_{-0.20}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f04647d55df3754efdfeeb54fc1e288b.svg?1327aa7b4c&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/> | <img alt="$15.28$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/829f235defecaee2215944ec4da79420.svg?6d55cfb700&invert_in_darkmode" align=middle width="37.3032pt" height="21.10812pt"/> |


#### Float <img alt="$N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?b2182879e4&invert_in_darkmode" align=middle width="24.21903pt" height="22.38192pt"/>, <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?ab43676d81&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/> and <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?4c5053aa82&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/>

![Disk only](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-z0-beta.png "Disk only")

#### Float <img alt="$N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?79f5c8e047&invert_in_darkmode" align=middle width="24.21903pt" height="22.38192pt"/>, <img alt="$N_B$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?1af8d2fd5&invert_in_darkmode" align=middle width="23.61183pt" height="22.38192pt"/>, <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?b3853a615b&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/> and <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?4be96edfa3&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/>

![Disk and bulge](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-nb-z0-beta.png "Disk and bulge")

#### Float <img alt="$N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?cf3bf446e8&invert_in_darkmode" align=middle width="24.21903pt" height="22.38192pt"/>, <img alt="$N_B$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?51cd7e0e60&invert_in_darkmode" align=middle width="23.61183pt" height="22.38192pt"/>, <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?5dd5123e8f&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/>, <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?cd6b1a9414&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/> and <img alt="$\alpha$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c745b9b57c145ec5577b82542b2df546.svg?316d1d0c52&invert_in_darkmode" align=middle width="10.537065pt" height="14.10255pt"/>

![Float alpha](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-nb-alpha-z0-beta.png "Float alpha")

 
### `data-corrected`

| <img alt="$N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?66dacbdc7f&invert_in_darkmode" align=middle width="24.21903pt" height="22.38192pt"/>             | <img alt="$N_B$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?b6061ebf0f&invert_in_darkmode" align=middle width="23.61183pt" height="22.38192pt"/>                     | <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?286b8f904d&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/>         | <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?3b80a16ca7&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/>       | <img alt="$\alpha$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c745b9b57c145ec5577b82542b2df546.svg?be5d836dae&invert_in_darkmode" align=middle width="10.537065pt" height="14.10255pt"/>      | TS |
|-------------------|---------------------------|---------------|---------------|---------------|----|
| <img alt="$1516977^{+559102}_{-445328}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/a2dbde4c286031fe29221e76b28d188b.svg?7c464e82be&invert_in_darkmode" align=middle width="106.72365pt" height="28.83969pt"/> | <img alt="$0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/29632a9bf827ce0200454dd32fc3be82.svg?a9132683de&invert_in_darkmode" align=middle width="8.188554pt" height="21.10812pt"/>                         | <img alt="$0.07^{+0.05}_{-0.03}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d205e9caf9a1f35c6a996c7a187b35e9.svg?da2f667a7d&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/> | <img alt="$2.09^{+0.07}_{-0.07}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/52ed15f3fcf5eaa6fb3287319b5599cc.svg?b3753ebe1b&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/>    | -             | <img alt="$0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/29632a9bf827ce0200454dd32fc3be82.svg?bcacf2aa19&invert_in_darkmode" align=middle width="8.188554pt" height="21.10812pt"/>  |
|<img alt="$1117715^{+484416}_{-363228}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f46a4868c7ddc11dd9f717ab7c69e06a.svg?dc9a184723&invert_in_darkmode" align=middle width="106.72365pt" height="28.83969pt"/>      | <img alt="$533365^{+529982}_{-246334}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8695d28150e98452914068400413fbf5.svg?631f1c6f7&invert_in_darkmode" align=middle width="98.53503pt" height="28.83969pt"/>    |  <img alt="$0.05^{+0.03}_{-0.02}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/31bc73fa691ca22a146f69461671e19e.svg?a63c1fed49&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/>| <img alt="$2.10^{+0.08}_{-0.07}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/9f056d70bba38553c351ba80f246f086.svg?63a2a97e93&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/>    | <img alt="$2.6$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/7509c44aea1098662d4f444907d5a427.svg?139e85077a&invert_in_darkmode" align=middle width="20.92629pt" height="21.10812pt"/>             | <img alt="$11.19$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/a92c8fcac4877df5c47adf945cb95ec6.svg?a09aef4caf&invert_in_darkmode" align=middle width="37.3032pt" height="21.10812pt"/> |
| <img alt="$1022473^{+463709}_{-331367}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/209b5e354269e9cf86774b47f449b347.svg?197f81fb7d&invert_in_darkmode" align=middle width="106.72365pt" height="28.83969pt"/> | <img alt="$1120457^{+1080567}_{-671964}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/03fa68f8631fec9cbd59c50469da095c.svg?df5b6f5e02&invert_in_darkmode" align=middle width="113.069715pt" height="28.83969pt"/> | <img alt="$0.06^{+0.04}_{-0.02}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/5eaf3c8b3988b60f7568e7a20f266d69.svg?11f79d9e5f&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/> | <img alt="$2.08^{+0.08}_{-0.07}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/32f07df8466c8593312dfe361b06cc27.svg?651831d91e&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/> | <img alt="$2.84^{+0.09}_{-0.25}$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/349ad4413b0561e4529c33ddd7186474.svg?142ed53979&invert_in_darkmode" align=middle width="62.82408pt" height="28.83969pt"/> | <img alt="$12.84$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/5b623fccabd46e08f79ddb507adddfe9.svg?890883914f&invert_in_darkmode" align=middle width="37.3032pt" height="21.10812pt"/> |


#### Float <img alt="$N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?1547e38252&invert_in_darkmode" align=middle width="24.21903pt" height="22.38192pt"/>, <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?4823d64685&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/> and <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?3a270c86a&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/>

![Disk only](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-z0-beta-data_eff.png "Disk only")

#### Float <img alt="$N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?db06e66aa&invert_in_darkmode" align=middle width="24.21903pt" height="22.38192pt"/>, <img alt="$N_B$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?c901c50ac5&invert_in_darkmode" align=middle width="23.61183pt" height="22.38192pt"/>, <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?e33441ca07&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/> and <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?e8ac02d511&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/>

![Disk and bulge](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-nb-z0-beta-data_eff.png "Disk and bulge")

#### Float <img alt="$N_D$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?aed0e1762a&invert_in_darkmode" align=middle width="24.21903pt" height="22.38192pt"/>, <img alt="$N_B$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?7f92527ce9&invert_in_darkmode" align=middle width="23.61183pt" height="22.38192pt"/>, <img alt="$z_0$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?c8cc6ed7&invert_in_darkmode" align=middle width="14.14413pt" height="14.10255pt"/>, <img alt="$\beta$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?7d8a642168&invert_in_darkmode" align=middle width="10.1277pt" height="22.74591pt"/> and <img alt="$\alpha$" src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c745b9b57c145ec5577b82542b2df546.svg?b9c38e1aa8&invert_in_darkmode" align=middle width="10.537065pt" height="14.10255pt"/>

![Float alpha](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/plots/2PC-nd-nb-alpha-z0-beta-data_eff.png "Float alpha")

 
