# GCE-2FIG
This code uses 2FIG sources to constrain spatial distributions of pulsar-like PS populations.  The main purpose of the code is to test the results presented in 1705.00009.

## Compiling and Running

This code package is written in `python` and `cython`.  To compile the `cython`, go into the `python/` subfolder and execute the file `make.sh`.  Note that you must have `cython` installed.  Examples of how to run the code are presented in the `examples/` subfolder.

## The likelihood Function

The `GCE-2FIG` code package uses the likelihood function presented in 1705.00009.  The likelihood is imlemented in the file `likelihood.pyx` in the `python/` subfolder.  The sky is spatially binned in a cartesian grid, and the data consists of the number of pulsar candidates detected in each bin.  The model parameters characterize the disk and Bulge point source populations, and as we scan over the model parameters we calculated the expected number of detected sources in each bin, utlizing the _Fermi_-provided efficiency function for detecting sources at a given flux and spatial position.  The likelihood is then given by the product over all pixels of the Poisson probabilities to observe the data counts given the predicted model counts.

As an option, we also incorporate the prior distribution in 1705.00009.

## The Expected Number PS Counts

Here, we describe how we calculate the expected number of PS counts, given the model parameters.  These computations are implemented in the file `get_counts_inline.pyx` in the `python/` subfolder.  First, we describe the computation for a general PS population, with density function <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e33a91659757dc5a5a317888c76bc940.svg?87a711d44a&invert_in_darkmode" align=middle width=57.31143pt height=24.56553pt/>, where <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e4c7375533c0164f8ecc75cfe4198ea0.svg?5e2024448e&invert_in_darkmode" align=middle width=48.84429pt height=24.56553pt/> are spatial coordinates with origin at the Earth, and luminosity function <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?37ce876f0d&invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/>.  Note that <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d30a65b936d8007addc9c789d5a7ae49.svg?e5a6b9bbb5&invert_in_darkmode" align=middle width=6.8238225pt height=22.74591pt/> and <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?7f2cead8e9&invert_in_darkmode" align=middle width=7.0284885pt height=22.74591pt/> are Galactic longitude and latitude, respectively, while <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?e40c7d7afe&invert_in_darkmode" align=middle width=7.6767405pt height=14.10255pt/> is the distance from the Earth.  The expected number of counts in a given spatial bin (labeled by <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4fe48dde86ac2d37419f0b35d57ac460.svg?7430ea8fe7&invert_in_darkmode" align=middle width=20.612625pt height=21.60213pt/>) and flux bin (labelled by <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/63bb9849783d01d91403bc9a5fea12a2.svg?c2b84c5ba3&invert_in_darkmode" align=middle width=9.041505pt height=22.74591pt/>) is given by the integral

<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/0b66f12996fe020e588077671462ab65.svg?6db75d50e1&invert_in_darkmode" align=middle width=384.76845pt height=47.505645pt/></p> 

where <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/2e71a536af409f25b2c3baebdba40859.svg?69543bf11f&invert_in_darkmode" align=middle width=35.22189pt height=14.10255pt/> is the efficiency factor for pular-like PS detection in that bin.  Note that we have assumed that <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/2e71a536af409f25b2c3baebdba40859.svg?9f76173c4c&invert_in_darkmode" align=middle width=35.22189pt height=14.10255pt/> is constant within the bin, so that we may bring it outside of the integral.
 

In the code, for both the disk and the Bulge, we normalize both <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?88e0b7f81a&invert_in_darkmode" align=middle width=8.46714pt height=14.10255pt/> and <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?e141113827&invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/> such that the full integral over all of space and flux is equal to the total number of sources <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f9c4988898e7f532b9f826a75014ed3c.svg?460833f7a9&invert_in_darkmode" align=middle width=14.94405pt height=22.38192pt/> for that population.  The number of sources is one of the model parameters.  To work with this convention, we must know how to properly normalize <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?be8d288134&invert_in_darkmode" align=middle width=8.46714pt height=14.10255pt/> and <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?ab3287125d&invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/> .

The normalization of <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?c87184bd09&invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/> is common to both the disk and Bulge.  Let 

<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/dc7e9cac24bdd8b9fb9aa0389e209134.svg?31778f06b&invert_in_darkmode" align=middle width=69.953235pt height=33.769395pt/></p>

for <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ddcb483302ed36a59286424aa5e0be17.svg?83c87ed8bf&invert_in_darkmode" align=middle width=11.14542pt height=22.38192pt/> in the range <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f8702ab460607cfc76e52498178edb88.svg?aedfcba9fe&invert_in_darkmode" align=middle width=86.361pt height=24.56553pt/> and be zero outside of this range.  We will normalize <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?b4d818c4a5&invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/> so that it integrates to unity, which implies that
 
<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/0a8a3c1792e9fdb071d6a6920dec93d7.svg?7b9f12204&invert_in_darkmode" align=middle width=147.489705pt height=41.283165pt/></p>  


We discuss <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?8a2fc36c1c&invert_in_darkmode" align=middle width=8.46714pt height=14.10255pt/> now for the disk and Bulge.  We will use the normalization that <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/91f9537316cdc8858af8eeb325936c77.svg?7f34dff96&invert_in_darkmode" align=middle width=84.237945pt height=26.70657pt/>.

### Bulge Sources 

We characterize the Bulge PS population by

<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/20a7cb8590592ee1ab7830a4d8541614.svg?dc41f9b4fc&invert_in_darkmode" align=middle width=90.27315pt height=33.5874pt/></p>

where <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/89f2e0d2d24bcf44db73aab8fc03252c.svg?e3475dd552&invert_in_darkmode" align=middle width=7.8435885pt height=14.10255pt/> is the distance from the Galactic Center, so long as <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/af61c60decc49acae8eeda9dca8fe898.svg?740cba9b51&invert_in_darkmode" align=middle width=55.277805pt height=17.65764pt/> .  <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?8ca56527d7&invert_in_darkmode" align=middle width=8.46714pt height=14.10255pt/> is equal to zero for larger values of <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/89f2e0d2d24bcf44db73aab8fc03252c.svg?d59c899f7c&invert_in_darkmode" align=middle width=7.8435885pt height=14.10255pt/> .  To have to proper normalization, it is straightforward to verify that
 
<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/9e79c22983590f5b42a2672ec8cf10a1.svg?dcb8fdc553&invert_in_darkmode" align=middle width=142.32603pt height=40.08807pt/></p>


### Disk Sources

The disk is modeled by the distribution 

<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ffd30829e06c1a1f6e2acec3a70ed61c.svg?56ab6e2bf3&invert_in_darkmode" align=middle width=192.0963pt height=18.569595pt/></p>

and we find that in this case

<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/a5bba1689089f869eba4d9fda47f21e7.svg?e22817b638&invert_in_darkmode" align=middle width=175.99395pt height=37.68171pt/></p>


### Code Implementation

Here, we describe in more detail how we implement the computations of the expected number of PS counts in the file `get_counts_inline.pyx`.

The cartesian grid is given by 12 linear-spaced bins in the inner 40<img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bda93e7eec1ea3bd03d7177c5b991481.svg?94abd8fa0f&invert_in_darkmode" align=middle width=6.7100715pt height=22.59873pt/> of the Milky Way:

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

return the number of counts in the Bulge in spatial <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d30a65b936d8007addc9c789d5a7ae49.svg?fa4da1999&invert_in_darkmode" align=middle width=6.8238225pt height=22.74591pt/> bin `i`, <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?17a1f238e2&invert_in_darkmode" align=middle width=7.0284885pt height=22.74591pt/> bin `j`, and flux bin `k`.  Most of the parameters above are defined previously.  However, there are a few parameters that are specific to the calculation method.  These include

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

Then, the bulk of the code proceeds through the loops used to perform the angular integration and the integration over <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?8a89a4c633&invert_in_darkmode" align=middle width=7.6767405pt height=14.10255pt/>:

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

Above, `i_ell` is the index over longitude, `i_b` is the index of latitude, and `l` is the index over the <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?a6455962f5&invert_in_darkmode" align=middle width=7.6767405pt height=14.10255pt/> bins.  In the first few lines we simply update the new values for `b` and `ell` and their cosines.  We then see if these values of `b` and `ell` are masked.  The statment `if incl > 0:` looks to see if we are too far away from the GC, such that we will completely miss the bulge.  This is relevant given the finite extent of the bulge, but this step is not needed for the disk contribution.  Then, we begin looping over the `s` integral, through the index `l`.  First, we find the true limits of the flux integral, given that fact that the luminosity function is cut off at `Lmin` and `Lmax`.  Then, we perform the integral over `L` to get `pref_L`.  The quantity `integral` contains both the integral over `s` and the integral over `L`.  

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

Unless floated, <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8d355ff26fea5f8b634222b72f19472f.svg?1d71b4dd35&invert_in_darkmode" align=middle width=60.814545pt height=21.10812pt/>,
<img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/5a341cc43b8429c286e732c54f035ae3.svg?579b791075&invert_in_darkmode" align=middle width=69.12081pt height=21.10812pt/>,
<img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/fa32ade5e29240b938ede1a9fe79cd17.svg?16a6101ecf&invert_in_darkmode" align=middle width=53.33328pt height=21.10812pt/>,
<img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/01a4d380c848b0b24c9d376c25ec7880.svg?bb6e383f3e&invert_in_darkmode" align=middle width=52.926885pt height=22.74591pt/>,
<img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f2ca6be66d1c9df73aac3a892fda5efe.svg?e247f3f01e&invert_in_darkmode" align=middle width=57.78663pt height=21.10812pt/>.

Floating parameters as in Table 2 of [1705.00009](https://arxiv.org/pdf/1705.00009.pdf):

| <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?7f98de2ea1&invert_in_darkmode" align=middle width=24.21903pt height=22.38192pt/>             | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?5c08000eeb&invert_in_darkmode" align=middle width=23.61183pt height=22.38192pt/>                     | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?390bab582e&invert_in_darkmode" align=middle width=14.14413pt height=14.10255pt/>         | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?95b13ea02f&invert_in_darkmode" align=middle width=10.1277pt height=22.74591pt/>       | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c745b9b57c145ec5577b82542b2df546.svg?e0fc56e9c5&invert_in_darkmode" align=middle width=10.537065pt height=14.10255pt/>      | TS |
|-------------------|---------------------------|---------------|---------------|---------------|----|
| <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ad0377333c6205a6c95d85ab794035a2.svg?e507dbddaf&invert_in_darkmode" align=middle width=98.53503pt height=28.83969pt/> | 0                         | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ea0db90a34d030368936947fffe939b6.svg?708131c72f&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/> | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6a9f8a50b291933d41f2a3053694952b.svg?379596a4e6&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/>    | -             | 0  |
| <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e163af18d2235279822f013c40633efc.svg?368b876e6e&invert_in_darkmode" align=middle width=92.00697pt height=28.83969pt/>         | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/3ac2d0ed96187aa38f78ba6f50b39bf4.svg?932eaa49f3&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/>| <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ba99f1a8f7eff41153feecabbfabf059.svg?274b92cef9&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/>    | -             | 6.5 |
| <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/68a76f1b147b020e4c248ce98fc00986.svg?ad58c86f6b&invert_in_darkmode" align=middle width=98.53503pt height=28.83969pt/> | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/2cc8900269e8130fc67df16485761672.svg?b7193d5ba8&invert_in_darkmode" align=middle width=98.53503pt height=28.83969pt/> | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/3ac2d0ed96187aa38f78ba6f50b39bf4.svg?37cd55d1c4&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/> | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ba99f1a8f7eff41153feecabbfabf059.svg?3fccb48db5&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/> | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/73e837ed10ca4f2bb41c14f1f114a1b9.svg?7576a783c1&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/> | 6.3 |







 
