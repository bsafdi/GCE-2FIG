###############################################################################
# get_counts.pyx
###############################################################################
#
# Return the predicted bulge counts in a bin (i,j,k), following the definition
# in 1705.00009
#
# Written: Nick Rodd, MIT, 5 August 2017
# Modified by Ben Safdi, UM August 15 2017
# 
###############################################################################

# Import basic functions
import numpy as np
cimport numpy as np
cimport cython

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

cdef double pi = np.pi
cdef double rodot = 8.5 # Earth GC distance in kpc
cdef double degtorad = pi/180.
cdef double fluxunits = 1.52529e37 # [converts MeV/cm^2 to erg/kpc^2]
# Below value is chosen for mix of speed and accuracy - not perfect
cdef double sbin = 0.05 # Width of integration bins in kpc

# These are the central bin values for l and b. Note l is positive to the left
cdef double angvals[12]
angvals = [18.33333333, 15., 11.66666667, 8.33333333, 5., 1.66666667, 
           -1.66666667, -5., -8.33333333, -11.66666667, -15., -18.33333333]

# These are the bin boundaries
cdef double[::1] ang_boundaries = np.linspace(-20.,20.,13,dtype=DTYPE)

# Bin edges for flux in units [MeV/cm^2/s]
cdef double fluxvals[9]
fluxvals = [1.00000000e-06, 1.46779927e-06, 2.15443469e-06, 3.16227766e-06, 
            4.64158883e-06, 6.81292069e-06, 1.00000000e-05, 3.16227766e-05, 
            1.00000000e-04]

cdef double area = pow(angvals[1]-angvals[0],2)*pow(2.*pi/360.,2.)


#########################################################
# This is New Code written by Ben Safdi                 # 
# Predicted Bulge PSRs in bin (i,j,k) [long, lat, flux] #
# Performs the angular integral instead of taking the   #
# value from the bin centre                             #
#########################################################
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef inline double Nbulge_full_ang_ijk(int i, int j, int k, double Nbulge, double omega_ijk, double alpha, double beta,double rcut , double Lmin, double Lmax ,int Ns ,int Nang, double theta_mask ) nogil:
    """ Return the number of bulge PSRs in bin (i,j,k) for given parameters
    ---arguments---
    Ns : number of points in s integration
    ----units---
    rcut : in kpc
    Lmin : erg s^{-1}
    Lmax : erg s^{-1}
    theta_mask : degrees, mask from the GC
    """
    # Determine angles and flux boundaries, convert flux to erg/kpc^2/s
    cdef double fluxmin = fluxvals[k] * fluxunits
    cdef double fluxmax = fluxvals[k+1] * fluxunits

    # Setup output and if variables
    cdef double Nijk = 0.
    cdef double prefactor, pref_rho, pref_L, lim_m, lim_p, integral, smin, smax,s,ds

    cdef int l

    cdef double incl 

    cdef double res = 0.0

    cdef int i_b, i_ell
    cdef double ell,b
    cdef double ell_start = ang_boundaries[i]
    cdef double b_start = ang_boundaries[j]
    cdef double d_ang = (ang_boundaries[1] - ang_boundaries[0])/float(Nang)

    cdef double total_res = 0.0

    cdef double coslval, cosbval

    ####
    ##cos theta for the mask
    cdef double costhetaval = cos(theta_mask * degtorad)


    # Common prefactors
    pref_rho = omega_ijk * Nbulge * (3. - alpha) / 4. / pi / pow(rcut,3 - alpha)
    pref_L = pow(4.*pi,1.-beta)*(pow(fluxmin, 1.-beta) - pow(fluxmax, 1.-beta))/(pow(Lmin,1-beta) - pow(Lmax,1-beta)) #(beta-1.)

    # for s integration
    smin = rodot - rcut
    smax = rodot + rcut
    ds = (smax - smin)/float(Ns)

    for i_b in range(Nang):
        b = b_start + i_b*d_ang
        if cosbval * coslval < costhetaval:
            for i_ell in range(Nang):
                ell = ell_start + i_ell*d_ang
                coslval = cos(ell * degtorad)
                cosbval = cos(b * degtorad)

                #res = 0.0
                incl = 1. - (1. - pow(rcut/rodot,2) ) * pow(coslval*cosbval,2)
                if incl > 0:
                    # pref_rho = omega_ijk * Nbulge * (3. - alpha) / 4. / pi / pow(rcut,3 - alpha)
                    # pref_L = pow(4.*pi,1.-beta)*(pow(fluxmin, 1.-beta) - pow(fluxmax, 1.-beta))/(pow(Lmin,1-beta) - pow(Lmax,1-beta)) #(beta-1.)

                    # smin = rodot - rcut
                    # smax = rodot + rcut
                    # ds = (smax - smin)/float(Ns)
                    integral = 0.0
                    s = smin
                    for l in range(Ns):
                        r_squared = (s*cosbval*coslval - rodot)**2 + s**2*(1. - pow(cosbval*coslval,2) )
                        #print r_squared
                        if r_squared < rcut**2:
                            integral += pow(s,4.-2.*beta)*pow(r_squared,-alpha/2.)
                        s += ds

                    total_res += cosbval * integral #* ds * pref_L * pref_rho
                    #total_res += res
            #print integral, ds, pref_L, pref_rho

    #print area, cosbval, integral, ds, pref_L, pref_rho
    return total_res * d_ang**2. * degtorad**2 * ds * pref_L * pref_rho



#########################################################
# This is New Code written by Ben Safdi                 # 
# Predicted Disk PSRs in bin (i,j,k) [long, lat, flux] #
# Performs the angular integral instead of taking the   #
# value from the bin centre                             #
#########################################################
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef inline double Ndisk_full_ang_ijk(int i, int j, int k, double Ndisk, double omega_ijk, double n, double sigma,double z0, double beta, double Lmin, double Lmax ,int Ns ,int Nang ,double smax, double theta_mask ) nogil:
    """ Return the number of bulge PSRs in bin (i,j,k) for given parameters
    ---arguments---
    Ns : number of points in s integration
    ----units---
    z0 : kpc
    sigma : kpc
    Lmin : erg s^{-1}
    Lmax : erg s^{-1}
    smax : 40 kpc default, how far to integrate out to from Sun?
    theta_mask : degrees, mask from the GC
    """
    # Determine angles and flux boundaries, convert flux to erg/kpc^2/s
    cdef double fluxmin = fluxvals[k] * fluxunits
    cdef double fluxmax = fluxvals[k+1] * fluxunits

    # Setup output and if variables
    cdef double Nijk = 0.
    cdef double prefactor, pref_rho, pref_L, lim_m, lim_p, integral,s,ds

    cdef int l

    cdef double res = 0.0

    cdef int i_b, i_ell
    cdef double ell,b
    cdef double ell_start = ang_boundaries[i]
    cdef double b_start = ang_boundaries[j]
    cdef double d_ang = (ang_boundaries[1] - ang_boundaries[0])/float(Nang)

    cdef double total_res = 0.0

    cdef double coslval, cosbval, sinbval

    ####
    ##cos theta for the mask
    cdef double costhetaval = cos(theta_mask * degtorad)

    cdef double R,z


    # Common prefactors
    pref_rho = omega_ijk * Ndisk /4./pi/z0/pow(sigma,n+2)/tgamma(n+2)
    pref_L = pow(4.*pi,1.-beta)*(pow(fluxmin, 1.-beta) - pow(fluxmax, 1.-beta))/(pow(Lmin,1-beta) - pow(Lmax,1-beta)) #(beta-1.)

    # for s integration
    smin = 0.0
    ds = (smax - smin)/float(Ns)

    for i_b in range(Nang):
        b = b_start + i_b*d_ang
        if cosbval * coslval < costhetaval:
            for i_ell in range(Nang):
                ell = ell_start + i_ell*d_ang
                coslval = cos(ell * degtorad)
                cosbval = cos(b * degtorad)
                #sinbval = sin(b * degtorad)

                res = 0.0
                integral = 0.0
                s = smin
                for l in range(Ns):
                    r_squared = (s*cosbval*coslval - rodot)**2 + s**2*(1. - pow(cosbval*coslval,2) )
                    z = s*sqrt(1 - cosbval**2)
                    R = sqrt( r_squared - z**2 )
                    #print r_squared
                    integral += pow(s,4.-2.*beta)*pow(R,n)*exp(-R/sigma - z/z0)
                    s += ds

                res = cosbval * integral #* ds * pref_L * pref_rho
                total_res += cosbval * integral
        #print integral, ds, pref_L, pref_rho

    #print area, cosbval, integral, ds, pref_L, pref_rho
    return total_res * d_ang**2. * degtorad**2 * ds * pref_L * pref_rho









