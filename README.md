# GCE-2FIG
This code uses 2FIG sources to constrain spatial distributions of pulsar-like PS populations.  The main purpose of the code is to test the results presented in 1705.00009.

## Compiling and Running

This code package is written in `python` and `cython`.  To compile the `cython`, go into the `python/` subfolder and execute the file `make.sh`.  Note that you must have `cython` installed.  Examples of how to run the code are presented in the `examples/` subfolder.

## The likelihood Function

The `GCE-2FIG` code package uses the likelihood function presented in 1705.00009.  The likelihood is imlemented in the file `likelihood.pyx` in the `python/` subfolder.  The sky is spatially binned in a cartesian grid, and the data consists of the number of pulsar candidates detected in each bin.  The model parameters characterize the disk and Bulge point source populations, and as we scan over the model parameters we calculated the expected number of detected sources in each bin, utlizing the _Fermi_-provided efficiency function for detecting sources at a given flux and spatial position.  The likelihood is then given by the product over all pixels of the Poisson probabilities to observe the data counts given the predicted model counts.

As an option, we also incorporate the prior distribution in 1705.00009.

## The Expected Number PS Counts

Here, we describe how we calculate the expected number of PS counts, given the model parameters.  These computations are implemented in the file `get_counts_inline.pyx` in the `python/` subfolder.  First, we describe the computation for a general PS population, with density function <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e33a91659757dc5a5a317888c76bc940.svg?invert_in_darkmode" align=middle width=57.31143pt height=24.56553pt/>, where <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/e4c7375533c0164f8ecc75cfe4198ea0.svg?invert_in_darkmode" align=middle width=48.84429pt height=24.56553pt/> are spatial coordinates with origin at the Earth, and luminosity function <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/>.  Note that <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/d30a65b936d8007addc9c789d5a7ae49.svg?invert_in_darkmode" align=middle width=6.8238225pt height=22.74591pt/> and <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?invert_in_darkmode" align=middle width=7.0284885pt height=22.74591pt/> are Galactic longitude and latitude, respectively, while <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?invert_in_darkmode" align=middle width=7.6767405pt height=14.10255pt/> is the distance from the Earth.  The expected number of counts in a given spatial bin (labeled by <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/4fe48dde86ac2d37419f0b35d57ac460.svg?invert_in_darkmode" align=middle width=20.612625pt height=21.60213pt/>) and flux bin (labelled by <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode" align=middle width=9.041505pt height=22.74591pt/>) is given by the integral
<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/0b66f12996fe020e588077671462ab65.svg?invert_in_darkmode" align=middle width=384.76845pt height=47.505645pt/></p> 
where <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/2e71a536af409f25b2c3baebdba40859.svg?invert_in_darkmode" align=middle width=35.22189pt height=14.10255pt/> is the efficiency factor for pular-like PS detection in that bin.  Note that we have assumed that <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/2e71a536af409f25b2c3baebdba40859.svg?invert_in_darkmode" align=middle width=35.22189pt height=14.10255pt/> is constant within the bin, so that we may bring it outside of the integral. 

In the code, for both the disk and the Bulge, we normalize both <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?invert_in_darkmode" align=middle width=8.46714pt height=14.10255pt/> and <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/> such that the full integral over all of space and flux is equal to the total number of sources <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode" align=middle width=14.94405pt height=22.38192pt/> for that population.  The number of sources is one of the model parameters.  To work with this convention, we must know how to properly normalize <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?invert_in_darkmode" align=middle width=8.46714pt height=14.10255pt/> and <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/>.

The normalization of <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/> is common to both the disk and Bulge.  Let 
<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/dc7e9cac24bdd8b9fb9aa0389e209134.svg?invert_in_darkmode" align=middle width=69.953235pt height=33.769395pt/></p>
for <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ddcb483302ed36a59286424aa5e0be17.svg?invert_in_darkmode" align=middle width=11.14542pt height=22.38192pt/> in the range <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/f8702ab460607cfc76e52498178edb88.svg?invert_in_darkmode" align=middle width=86.361pt height=24.56553pt/> and be zero outside of this range.  We will normalize <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bb7f89046aaec2638cc892ac2d6b7b12.svg?invert_in_darkmode" align=middle width=49.96266pt height=24.56553pt/> so that it integrates to unity, which implies that 
<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/0a8a3c1792e9fdb071d6a6920dec93d7.svg?invert_in_darkmode" align=middle width=147.489705pt height=41.283165pt/></p>  

We discuss <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?invert_in_darkmode" align=middle width=8.46714pt height=14.10255pt/> now for the disk and Bulge.  We will use the normalization that <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/91f9537316cdc8858af8eeb325936c77.svg?invert_in_darkmode" align=middle width=84.237945pt height=26.70657pt/>.

### Bulge Sources 

We characterize the Bulge PS population by
<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/20a7cb8590592ee1ab7830a4d8541614.svg?invert_in_darkmode" align=middle width=90.27315pt height=33.5874pt/></p>
where <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/89f2e0d2d24bcf44db73aab8fc03252c.svg?invert_in_darkmode" align=middle width=7.8435885pt height=14.10255pt/> is the distance from the Galactic Center, so long as <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/af61c60decc49acae8eeda9dca8fe898.svg?invert_in_darkmode" align=middle width=55.277805pt height=17.65764pt/>.  <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?invert_in_darkmode" align=middle width=8.46714pt height=14.10255pt/> is equal to zero for larger values of <img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/89f2e0d2d24bcf44db73aab8fc03252c.svg?invert_in_darkmode" align=middle width=7.8435885pt height=14.10255pt/>.  To have to proper normalization, it is straightforward to verify that 
<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/9e79c22983590f5b42a2672ec8cf10a1.svg?invert_in_darkmode" align=middle width=142.32603pt height=40.08807pt/></p>

### Disk Sources

The disk is modeled by the distribution 
<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/ffd30829e06c1a1f6e2acec3a70ed61c.svg?invert_in_darkmode" align=middle width=192.0963pt height=18.569595pt/></p>
and we find that in this case
<p align="center"><img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/a5bba1689089f869eba4d9fda47f21e7.svg?invert_in_darkmode" align=middle width=175.99395pt height=37.68171pt/></p>

### Code Implementation

Here, we describe in more detail how we implement the computations of the expected number of PS counts in the file `get_counts_inline.pyx`.

The cartesian grid is given by 12 linear-spaced bins in the inner 40<img src="https://rawgit.com/bsafdi/GCE-2FIG/master/svgs/bda93e7eec1ea3bd03d7177c5b991481.svg?invert_in_darkmode" align=middle width=6.7100715pt height=22.59873pt/> of the Milky Way:
```python
cdef double[::1] ang_boundaries = np.linspace(-20.,20.,13,dtype=DTYPE)
``` 





 
