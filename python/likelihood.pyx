###############################################################################
# likelihood.pyx
###############################################################################
#
# Evaluate the likelihood as defined in 1705.00009.
#
#
# Written: Siddharth Mishra-Sharma, Nick Rodd, and Ben Safdi 15 August 2017
#
###############################################################################


# Import basic functions
import numpy as np
cimport numpy as np
cimport cython


# import get_counts as gc
cimport get_counts_inline as gc

DTYPE = np.float
ctypedef np.float_t DTYPE_t

# C functions
cdef extern from "math.h":
    double pow(double x, double y) nogil
    double sqrt(double x) nogil
    double cos(double x) nogil
    double tgamma(double x) nogil
    double exp(double x) nogil
    double sin(double x) nogil
    double fabs(double x) nogil
    double tan(double x) nogil

    double log(double x) nogil
    double lgamma(double x) nogil

cdef double pi = np.pi



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef double ll(double[:, :, ::1] PSR_data, double[:, :, ::1] omega_ijk, 
                double Nbulge, double Ndisk, double alpha, double n, 
                double sigma, double z0, double Lmax_disk, double Lmax_bulge, double beta_disk, double beta_bulge, double rcut = 3.0, double Lmin = 1.0e31,int Ns = 200, int Nang = 1,double smax_disk = 30,double theta_mask=0.0,
                int use_prior = 0, int Nang_prior = 20) nogil:
    """ Calculate the likelihood as a function of the bulge and disk params
    theta_mask in degrees 
    """

    # Setup loop variables
    cdef double omega_val, Nmodel, Nobs, N_bulge, N_disk
    cdef Py_ssize_t i, j, k
    
    # Now calculate the likelihood over all bins
    cdef double ll = 0.
    for i in range(12): # loop over longitude
        for j in range(12): # loop over latitude
            for k in range(8): # loop over flux
                if i != 6:
                    omega_val = omega_ijk[i,j,k] 

                    N_disk = gc.Ndisk_full_ang_ijk(i,j,k,Ndisk, omega_val,n,sigma,z0,beta_disk,Lmin,Lmax_disk,Ns,Nang,smax_disk,theta_mask)
                    N_bulge = gc.Nbulge_full_ang_ijk(i,j,k,Nbulge, omega_val,alpha,beta_bulge,rcut, Lmin,Lmax_bulge,Ns,Nang,theta_mask)

                    Nmodel = N_disk + N_bulge

                    Nobs = PSR_data[i,j,k]


                    ll += Nobs*log(Nmodel) - Nmodel #- lgamma(Nobs + 1)
             
    cdef double Nmodeltot = 0.
    cdef double pmid = 174.
    cdef double psig = 63.

    if use_prior:  # If using a prior on total number of sources

        # Add bulge
        Nmodeltot += gc.Nbulge_total(Nbulge,  alpha,  beta_bulge, rcut ,  Lmin,  Lmax_bulge , Ns , Nang_prior,  theta_mask)

        # Add disk
        Nmodeltot += gc.Ndisk_total(Ndisk,  n,  sigma, z0,  beta_disk,  Lmin,  Lmax_disk , Ns , Nang_prior , smax_disk,  theta_mask )
        
        # Now add the prior term - note I'm sure they have a sign error in that too
        ll -= pow(Nmodeltot - pmid, 2.) / (2. * pow(psig, 2.))

    return ll



