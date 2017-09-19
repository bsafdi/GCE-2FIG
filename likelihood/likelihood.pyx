###############################################################################
# likelihood.pyx
###############################################################################
#
# Evaluate the likelihood as defined in 1705.00009, using our custom PSR count
# functions in PSR_counts.pyx 
#
###############################################################################


# Import functions including custom PSR count functions
import numpy as np
cimport numpy as np
cimport cython
cimport PSR_counts as pc

# C functions
cdef extern from "math.h":
    double pow(double x, double y) nogil
    double log(double x) nogil


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef double ll(double[:, :, ::1] PSR_data, double[:, :, ::1] omega_ijk, 
                double Nbulge, double Ndisk, double alpha, double n, 
                double sigma, double z0, double Lmax_disk, double Lmax_bulge, 
                double beta_disk, double beta_bulge, double rcut = 3., 
                double Lmin = 1.e31, int Ns = 200, int Nang = 1, 
                double smax_disk = 30, double theta_mask = 2.,
                int use_prior = 0, int Nang_prior = 20) nogil:
    """ Calculate the likelihood as a function of the bulge and disk params

    ---arguments---
    PSR_data : data binned in [longitude, latitude, flux]
    omega_ijk : efficiency in same binning as data
    Nbulge : number of bulge PSRs
    Ndisk : number of disk PSRs
    alpha : inner slope of the bulge spatial profile
    n : cylindrical radial index of the disk profile
    sigma : exponential cutoff of disk radius [kpc]
    z0 : exponential cutoff of disk height [kpc]
    Lmax_disk : max luminosity of the disk [erg/s]
    Lmax_bulge : max luminosity of the bulge [erg/s]
    beta_disk : index for luminosity function of disk
    beta_bulge : index for luminosity function of bulge
    Nang : number of angular steps per bin for the likelihood integral
    smax_disk : maximum distance to integrate through the disk from Earth [kpc]
    theta_mask : angle from the GC to be masked
    use_prior : whether to use the prior
    Nang_prior : number of angular steps per bin for the prior integral
    """

    # Setup loop variables
    cdef double Nobs, omega_val, Nmodel
    cdef Py_ssize_t i, j, k
    
    # Now calculate the log likelihood by summing over all bins
    cdef double ll = 0.
    # Data is binned from top left (l,b)=(+20,+20) to bottom right (-20,-20)
    for i in range(12): # loop over longitude
        for j in range(12): # loop over latitude
            for k in range(8): # loop over flux
                Nobs = PSR_data[i,j,k]
                omega_val = omega_ijk[i,j,k] 

                Nmodel = 0.
                Nmodel += pc.Ndisk_ijk(i,j,k,Ndisk,omega_val,n,sigma,z0,
                                      beta_disk,Lmin,Lmax_disk,Ns,Nang,
                                      smax_disk,theta_mask)
                Nmodel += pc.Nbulge_ijk(i,j,k,Nbulge,omega_val,alpha,
                                        beta_bulge,rcut, Lmin,Lmax_bulge,Ns,
                                        Nang,theta_mask)

                ll += Nobs*log(Nmodel) - Nmodel 
             
    # Add the prior contribution if specified
    cdef double Nmodeltot = 0.
    cdef double pmid = 174.
    cdef double psig = 63.

    if use_prior:
        # By default we have no mask in the prior, so hard coded to 0
        Nmodeltot += pc.Nbulge_total(Nbulge,alpha,beta_bulge,rcut,Lmin,
                                     Lmax_bulge,Ns,Nang_prior,0.)
        Nmodeltot += pc.Ndisk_total(Ndisk,n,sigma,z0,beta_disk,Lmin,Lmax_disk,
                                    Ns,Nang_prior,smax_disk,0.)
        
        # Add the prior contribution to the log likelihood
        ll -= pow(Nmodeltot - pmid, 2.) / (2. * pow(psig, 2.))

    return ll
