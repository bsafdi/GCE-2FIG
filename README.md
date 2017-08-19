# GCE-2FIG
This code uses 2FIG sources to constrain spatial distributions of pulsar-like PS populations.  The main purpose of the code is to test the results presented in 1705.00009.

## Compiling and Running

This code package is written in `python` and `cython`.  To compile the `cython`, go into the `python/` subfolder and execute the file `make.sh`.  Note that you must have `cython` installed.  Examples of how to run the code are presented in the `examples/` subfolder.

## The likelihood Function

The `GCE-2FIG` code package uses the likelihood function presented in 1705.00009.  The likelihood is imlemented in the file `likelihood.pyx` in the `python/` subfolder.  The sky is spatially binned in a cartesian grid, and the data consists of the number of pulsar candidates detected in each bin.  The model parameters characterize the disk and Bulge point source populations, and as we scan over the model parameters we calculated the expected number of detected sources in each bin, utlizing the _Fermi_-provided efficiency function for detecting sources at a given flux and spatial position.  The likelihood is then given by the product over all pixels of the Poisson probabilities to observe the data counts given the predicted model counts.

As an option, we also incorporate the prior distribution in 1705.00009.

## The Expected Number PS Counts

Here, we describe how we calculate the expected number of PS counts, given the model parameters.  These computations are implemented in the file `get_counts_inline.pyx` in the `python/` subfolder.  First, we describe the computation for a general PS population, with density function <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e33a91659757dc5a5a317888c76bc940.svg?ab22ebd05&invert_in_darkmode" align=middle width=57.31143pt height=24.56553pt/>, where <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e4c7375533c0164f8ecc75cfe4198ea0.svg?285cec6fd&invert_in_darkmode" align=middle width=48.84429pt height=24.56553pt/> are spatial coordinates with origin at the Earth, and luminosity function <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?3c06e70f3a&invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/>.  Note that <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d30a65b936d8007addc9c789d5a7ae49.svg?3dcf459d4&invert_in_darkmode" align=middle width=6.8238225pt height=22.74591pt/> and <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?c68d6d5887&invert_in_darkmode" align=middle width=7.0284885pt height=22.74591pt/> are Galactic longitude and latitude, respectively, while <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?3675391101&invert_in_darkmode" align=middle width=7.6767405pt height=14.10255pt/> is the distance from the Earth.  The expected number of counts in a given spatial bin (labeled by <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4fe48dde86ac2d37419f0b35d57ac460.svg?2acecebf12&invert_in_darkmode" align=middle width=20.612625pt height=21.60213pt/>) and flux bin (labelled by <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/63bb9849783d01d91403bc9a5fea12a2.svg?be96b1668&invert_in_darkmode" align=middle width=9.041505pt height=22.74591pt/>) is given by the integral

<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/0b66f12996fe020e588077671462ab65.svg?9bb505bebc&invert_in_darkmode" align=middle width=384.76845pt height=47.505645pt/></p> 

where <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/2e71a536af409f25b2c3baebdba40859.svg?1cb9a216a&invert_in_darkmode" align=middle width=35.22189pt height=14.10255pt/> is the efficiency factor for pular-like PS detection in that bin.  Note that we have assumed that <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/2e71a536af409f25b2c3baebdba40859.svg?23b3e98049&invert_in_darkmode" align=middle width=35.22189pt height=14.10255pt/> is constant within the bin, so that we may bring it outside of the integral.
 

In the code, for both the disk and the Bulge, we normalize both <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?8124a88f40&invert_in_darkmode" align=middle width=8.46714pt height=14.10255pt/> and <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?7b891cb6dd&invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/> such that the full integral over all of space and flux is equal to the total number of sources <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f9c4988898e7f532b9f826a75014ed3c.svg?dbc9c45733&invert_in_darkmode" align=middle width=14.94405pt height=22.38192pt/> for that population.  The number of sources is one of the model parameters.  To work with this convention, we must know how to properly normalize <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?8530b0d3b1&invert_in_darkmode" align=middle width=8.46714pt height=14.10255pt/> and <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?dd1d8ced9f&invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/> .

The normalization of <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?7f843fd378&invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/> is common to both the disk and Bulge.  Let 

<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/dc7e9cac24bdd8b9fb9aa0389e209134.svg?a1e6c1926d&invert_in_darkmode" align=middle width=69.953235pt height=33.769395pt/></p>

for <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ddcb483302ed36a59286424aa5e0be17.svg?c66b413f7f&invert_in_darkmode" align=middle width=11.14542pt height=22.38192pt/> in the range <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f8702ab460607cfc76e52498178edb88.svg?e07c74a8fc&invert_in_darkmode" align=middle width=86.361pt height=24.56553pt/> and be zero outside of this range.  We will normalize <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?5ba9a9c76a&invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/> so that it integrates to unity, which implies that
 
<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/0a8a3c1792e9fdb071d6a6920dec93d7.svg?746ab08ce4&invert_in_darkmode" align=middle width=147.489705pt height=41.283165pt/></p>  


We discuss <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?de87b0790d&invert_in_darkmode" align=middle width=8.46714pt height=14.10255pt/> now for the disk and Bulge.  We will use the normalization that <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/91f9537316cdc8858af8eeb325936c77.svg?127444f38a&invert_in_darkmode" align=middle width=84.237945pt height=26.70657pt/>.

### Bulge Sources 

We characterize the Bulge PS population by

<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/20a7cb8590592ee1ab7830a4d8541614.svg?a311de0c3d&invert_in_darkmode" align=middle width=90.27315pt height=33.5874pt/></p>

where <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/89f2e0d2d24bcf44db73aab8fc03252c.svg?27c088552f&invert_in_darkmode" align=middle width=7.8435885pt height=14.10255pt/> is the distance from the Galactic Center, so long as <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/af61c60decc49acae8eeda9dca8fe898.svg?9122486d5&invert_in_darkmode" align=middle width=55.277805pt height=17.65764pt/> .  <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?a4b4891757&invert_in_darkmode" align=middle width=8.46714pt height=14.10255pt/> is equal to zero for larger values of <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/89f2e0d2d24bcf44db73aab8fc03252c.svg?5c07290d9b&invert_in_darkmode" align=middle width=7.8435885pt height=14.10255pt/> .  To have to proper normalization, it is straightforward to verify that
 
<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/9e79c22983590f5b42a2672ec8cf10a1.svg?ae44e7caa0&invert_in_darkmode" align=middle width=142.32603pt height=40.08807pt/></p>


### Disk Sources

The disk is modeled by the distribution 

<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ffd30829e06c1a1f6e2acec3a70ed61c.svg?60b1cf013a&invert_in_darkmode" align=middle width=192.0963pt height=18.569595pt/></p>

and we find that in this case

<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/a5bba1689089f869eba4d9fda47f21e7.svg?3afe1b633&invert_in_darkmode" align=middle width=175.99395pt height=37.68171pt/></p>


### Code Implementation

Here, we describe in more detail how we implement the computations of the expected number of PS counts in the file `get_counts_inline.pyx`.

The cartesian grid is given by 12 linear-spaced bins in the inner 40<img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bda93e7eec1ea3bd03d7177c5b991481.svg?532e37527d&invert_in_darkmode" align=middle width=6.7100715pt height=22.59873pt/> of the Milky Way:

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

return the number of counts in the Bulge in spatial <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d30a65b936d8007addc9c789d5a7ae49.svg?3c2f6ebb54&invert_in_darkmode" align=middle width=6.8238225pt height=22.74591pt/> bin `i`, <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?7d5a3d4bba&invert_in_darkmode" align=middle width=7.0284885pt height=22.74591pt/> bin `j`, and flux bin `k`.  Most of the parameters above are defined previously.  However, there are a few parameters that are specific to the calculation method.  These include

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

Then, the bulk of the code proceeds through the loops used to perform the angular integration and the integration over <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?4b4ef10587&invert_in_darkmode" align=middle width=7.6767405pt height=14.10255pt/>:

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

Above, `i_ell` is the index over longitude, `i_b` is the index of latitude, and `l` is the index over the <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?2ad2861f9d&invert_in_darkmode" align=middle width=7.6767405pt height=14.10255pt/> bins.  In the first few lines we simply update the new values for `b` and `ell` and their cosines.  We then see if these values of `b` and `ell` are masked.  The statment `if incl > 0:` looks to see if we are too far away from the GC, such that we will completely miss the bulge.  This is relevant given the finite extent of the bulge, but this step is not needed for the disk contribution.  Then, we begin looping over the `s` integral, through the index `l`.  First, we find the true limits of the flux integral, given that fact that the luminosity function is cut off at `Lmin` and `Lmax`.  Then, we perform the integral over `L` to get `pref_L`.  The quantity `integral` contains both the integral over `s` and the integral over `L`.  

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

Fits performed with `MultiNest`. 

Unless floated, values assumed:
- <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8d355ff26fea5f8b634222b72f19472f.svg?a87d82957b&invert_in_darkmode" align=middle width=60.814545pt height=21.10812pt/>,
- <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/5a341cc43b8429c286e732c54f035ae3.svg?7798c49311&invert_in_darkmode" align=middle width=69.12081pt height=21.10812pt/>,
- <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/fa32ade5e29240b938ede1a9fe79cd17.svg?c1f11312ec&invert_in_darkmode" align=middle width=53.33328pt height=21.10812pt/>,
- <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/01a4d380c848b0b24c9d376c25ec7880.svg?d65fdef281&invert_in_darkmode" align=middle width=52.926885pt height=22.74591pt/>,
- <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f2ca6be66d1c9df73aac3a892fda5efe.svg?1a5515053&invert_in_darkmode" align=middle width=57.78663pt height=21.10812pt/>.

Priors used: 
- <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/baca37ec0463c58e9b4bca0c53e033d2.svg?ba5447ed56&invert_in_darkmode" align=middle width=55.98087pt height=22.38192pt/>: [0,3000000]
- <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c745b9b57c145ec5577b82542b2df546.svg?c42056cc6e&invert_in_darkmode" align=middle width=10.537065pt height=14.10255pt/>: [2.1, 5.0]
- <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?c6bc35415f&invert_in_darkmode" align=middle width=14.14413pt height=14.10255pt/>: [0.01, 2.0]
- <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?d52c094f12&invert_in_darkmode" align=middle width=10.1277pt height=22.74591pt/>: [1.1,3.0]

Floating parameters as in Table 2 of [1705.00009](https://arxiv.org/pdf/1705.00009.pdf):

| <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?7860a0e6e1&invert_in_darkmode" align=middle width=24.21903pt height=22.38192pt/>             | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?9ba7e058d0&invert_in_darkmode" align=middle width=23.61183pt height=22.38192pt/>                     | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?5ee11657d9&invert_in_darkmode" align=middle width=14.14413pt height=14.10255pt/>         | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?6b7bcbf1f4&invert_in_darkmode" align=middle width=10.1277pt height=22.74591pt/>       | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c745b9b57c145ec5577b82542b2df546.svg?17f527edff&invert_in_darkmode" align=middle width=10.537065pt height=14.10255pt/>      | TS |
|-------------------|---------------------------|---------------|---------------|---------------|----|
| <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/a3222eb0b8a032f09b5ca68e5e2d71b8.svg?65ebbe7926&invert_in_darkmode" align=middle width=106.72365pt height=28.83969pt/> | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/29632a9bf827ce0200454dd32fc3be82.svg?a579594497&invert_in_darkmode" align=middle width=8.188554pt height=21.10812pt/>                         | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d205e9caf9a1f35c6a996c7a187b35e9.svg?30b3d761a0&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/> | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/127ad06b2f283a2f150297d49a237862.svg?7199efe2cb&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/>    | -             | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/29632a9bf827ce0200454dd32fc3be82.svg?2889a731a1&invert_in_darkmode" align=middle width=8.188554pt height=21.10812pt/>  |
|<img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/040e105630eeaca3afefb963ce5a5984.svg?d54ebdc743&invert_in_darkmode" align=middle width=106.72365pt height=28.83969pt/>      | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb10ef89e55b4b8fa086bcc19284c79c.svg?43c310a65b&invert_in_darkmode" align=middle width=98.53503pt height=28.83969pt/>    |  <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/31bc73fa691ca22a146f69461671e19e.svg?1ff96ff87&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/>| <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d90b994787399cb633c9787930ea00d6.svg?8c94c5c037&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/>    | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/7509c44aea1098662d4f444907d5a427.svg?1f64d99f93&invert_in_darkmode" align=middle width=20.92629pt height=21.10812pt/>             | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/edce064f9c22a4cc7aa63ea53e737a63.svg?5ea55d9b50&invert_in_darkmode" align=middle width=29.114745pt height=21.10812pt/> |
| <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/7b6cf9ff1431d0616b9e58f8b6c524bc.svg?6eedc229d2&invert_in_darkmode" align=middle width=106.72365pt height=28.83969pt/> | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/fc625dd50232b817c2bb316bbba41859.svg?af476ed79f&invert_in_darkmode" align=middle width=104.88126pt height=28.83969pt/> | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/5eaf3c8b3988b60f7568e7a20f266d69.svg?e37a9e224f&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/> | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c706cf7a9ce2586d0d47aae399e4202c.svg?e463f2717c&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/> | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bbe989c4b1e7c5e84148d56a17cf16bc.svg?d352da269e&invert_in_darkmode" align=middle width=62.82408pt height=28.83969pt/> | <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4c6ff5d6a1bbf365a109b7c3b6a3a148.svg?73b437fcba&invert_in_darkmode" align=middle width=29.114745pt height=21.10812pt/> |


### Float <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?13ca03d6d0&invert_in_darkmode" align=middle width=24.21903pt height=22.38192pt/>, <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?d3282ffbcd&invert_in_darkmode" align=middle width=14.14413pt height=14.10255pt/> and <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?54b0d72cc8&invert_in_darkmode" align=middle width=10.1277pt height=22.74591pt/>

![Disk only](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/nd.png "Disk only")

### Float <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?9b3fa43fe8&invert_in_darkmode" align=middle width=24.21903pt height=22.38192pt/>, <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?3f9fdae523&invert_in_darkmode" align=middle width=23.61183pt height=22.38192pt/>, <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?36f4f7a365&invert_in_darkmode" align=middle width=14.14413pt height=14.10255pt/> and <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?e70f79cb5e&invert_in_darkmode" align=middle width=10.1277pt height=22.74591pt/>

![Disk and bulge](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/ndnb.png "Disk and bulge")

### Float <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e1337360d9fc04449e0bd6a786b5ad20.svg?1290e84ebf&invert_in_darkmode" align=middle width=24.21903pt height=22.38192pt/>, <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f7a0c7db2f0eff424a274f52f4a5c8d6.svg?c86e88c52c&invert_in_darkmode" align=middle width=23.61183pt height=22.38192pt/>, <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d1a81d9dc6dd30e43ba27c5490a34a32.svg?a8f06b784a&invert_in_darkmode" align=middle width=14.14413pt height=14.10255pt/>, <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?ba938c3bb7&invert_in_darkmode" align=middle width=10.1277pt height=22.74591pt/> and <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/c745b9b57c145ec5577b82542b2df546.svg?77dec16e23&invert_in_darkmode" align=middle width=10.537065pt height=14.10255pt/>

![Float alpha](https://raw.githubusercontent.com/bsafdi/GCE-2FIG/master/examples/ndnbalpha.png "Float alpha")



 
