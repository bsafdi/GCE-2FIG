###############################################################################
# get_counts_inline.pyx
###############################################################################
#
# Return the predicted bulge and disk counts in a bin (i,j,k), 
# following the definition in 1705.00009.  Also compute the number of counts 
# above in a given flux range over the full sky for use with the prior
# these functions are called by likelihood.pyx
#
# Written: Siddharth Mishra-Sharma, Nick Rodd, and Ben Safdi 15 August 2017
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

# Useful numbres
cdef double pi = np.pi
cdef double rodot = 8.5 # Earth GC distance in kpc
cdef double degtorad = pi/180.
cdef double fluxunits = 1.52529e37 # [converts MeV/cm^2 to erg/kpc^2]
# Below value is chosen for mix of speed and accuracy - not perfect
#cdef double sbin = 0.05 # Width of integration bins in kpc

# # These are the central bin values for l and b. Note l is positive to the left
# cdef double angvals[12]
# angvals = [18.33333333, 15., 11.66666667, 8.33333333, 5., 1.66666667, 
#            -1.66666667, -5., -8.33333333, -11.66666667, -15., -18.33333333]

# These are the angular bin boundaries in degrees for l and b
cdef double[::1] ang_boundaries = np.linspace(-20.,20.,13,dtype=DTYPE)

# Bin edges for flux in units [MeV/cm^2/s]
cdef double fluxvals[9]
fluxvals = [1.00000000e-06, 1.46779927e-06, 2.15443469e-06, 3.16227766e-06, 
            4.64158883e-06, 6.81292069e-06, 1.00000000e-05, 3.16227766e-05, 
            1.00000000e-04]

# # Pixel area.
# cdef double area = pow(angvals[1]-angvals[0],2)*pow(2.*pi/360.,2.)


#########################################################
# Predicted Bulge PSRs in bin (i,j,k) [long, lat, flux] #
# Performs the angular integral instead of taking the   #
# value from the bin center                             #
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
    # pref_L = pow(4.*pi,1.-beta)*(pow(fluxmin, 1.-beta) - pow(fluxmax, 1.-beta))/(pow(Lmin,1-beta) - pow(Lmax,1-beta)) #(beta-1.)
    cdef double pref_L_norm = 1./(pow(Lmin,1-beta) - pow(Lmax,1-beta))

    cdef double flux_min_eval, flux_max_eval

    # for s integration
    smin = rodot - rcut
    smax = rodot + rcut
    ds = (smax - smin)/float(Ns)

    # loop over angles for angular integral
    for i_b in range(Nang):
        b = b_start + i_b*d_ang + d_ang/2.
        for i_ell in range(Nang):
            ell = ell_start + i_ell*d_ang + d_ang/2.
            coslval = cos(ell * degtorad)
            cosbval = cos(b * degtorad)

            if cosbval * coslval < costhetaval: #is it masked?

                incl = 1. - (1. - pow(rcut/rodot,2) ) * pow(coslval*cosbval,2)
                if incl > 0: # should be bother, or is the answer going to be zero?
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

    return total_res * d_ang**2. * degtorad**2 * ds * pref_rho



######################################################### 
# Predicted Disk PSRs in bin (i,j,k) [long, lat, flux]  #
# Performs the angular integral instead of taking the   #
# value from the bin center                             #
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
    #pref_L = pow(4.*pi,1.-beta)*(pow(fluxmin, 1.-beta) - pow(fluxmax, 1.-beta))/(pow(Lmin,1-beta) - pow(Lmax,1-beta)) 
    cdef double pref_L_norm = 1./(pow(Lmin,1-beta) - pow(Lmax,1-beta))

    cdef double flux_min_eval, flux_max_eval

    # for s integration
    smin = 0.01
    ds = (smax - smin)/float(Ns)

    # loop over angles
    for i_b in range(Nang):
        b = b_start + i_b*d_ang + d_ang/2.
        for i_ell in range(Nang):
            ell = ell_start + i_ell*d_ang + d_ang/2.
            coslval = cos(ell * degtorad)
            cosbval = cos(b * degtorad)
            if cosbval * coslval < costhetaval: #is it masked?

                res = 0.0
                integral = 0.0
                s = smin
                for l in range(Ns):

                    flux_min_eval = max(Lmin/s**2,4*pi*fluxmin)
                    flux_max_eval = min(Lmax/s**2,4*pi*fluxmax)

                    if flux_max_eval < flux_min_eval:
                        pref_L = 0.0
                    else:
                        pref_L = (pow(flux_min_eval, 1.-beta) - pow(flux_max_eval, 1.-beta))*pref_L_norm

                    r_squared = (s*cosbval*coslval - rodot)**2 + s**2*(1. - pow(cosbval*coslval,2) )
                    z = s*sqrt(1 - cosbval**2)
                    R = sqrt( r_squared - z**2 )
                    integral += pref_L * pow(s,4.-2.*beta)*pow(R,n)*exp(-R/sigma - z/z0)
                    s += ds

                res = cosbval * integral 
                total_res += cosbval * integral

    return total_res * d_ang**2. * degtorad**2 * ds  * pref_rho


#########################################################
#########################################################
#########################################################
# Below are the functions that calculate the number of  #
# counts for the prior                                  #
#########################################################
#########################################################
#########################################################



#########################################################
# Predicted Bulge PSRs, integraded over angs            # 
#########################################################
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef inline double Nbulge_total(double Nbulge, double alpha, double beta,double rcut , double Lmin, double Lmax ,int Ns ,int Nang, double theta_mask) nogil:
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
    cdef double fluxmin = 1.8e-5 * fluxunits #MeV cm^{-2} s^{-1}
    cdef double fluxmax = 0.003896 * fluxunits #MeV cm^{-2} s^{-1}

    # Setup output and if variables
    cdef double Nijk = 0.
    cdef double prefactor, pref_rho, pref_L, lim_m, lim_p, integral, smin, smax,s,ds

    cdef int l

    cdef double incl 

    cdef double res = 0.0

    cdef int i_b, i_ell
    cdef double ell,b
    cdef double ell_start = -15.0 # we will integrate for -15 to +15 in lat and lon
    cdef double b_start = -15.0 
    cdef double d_ang = (30.0)/float(Nang)

    cdef double total_res = 0.0

    cdef double coslval, cosbval

    ####
    ##cos theta for the mask
    cdef double costhetaval = cos(theta_mask * degtorad)


    # Common prefactors
    pref_rho = Nbulge * (3. - alpha) / 4. / pi / pow(rcut,3 - alpha)
    # pref_L = pow(4.*pi,1.-beta)*(pow(fluxmin, 1.-beta) - pow(fluxmax, 1.-beta))/(pow(Lmin,1-beta) - pow(Lmax,1-beta)) #(beta-1.)
    cdef double pref_L_norm = 1./(pow(Lmin,1-beta) - pow(Lmax,1-beta))

    cdef double flux_min_eval, flux_max_eval

    # for s integration
    smin = rodot - rcut
    smax = rodot + rcut
    ds = (smax - smin)/float(Ns)

    for i_b in range(Nang):
        b = b_start + i_b*d_ang + d_ang/2.
        for i_ell in range(Nang):
            ell = ell_start + i_ell*d_ang + d_ang/2.
            coslval = cos(ell * degtorad)
            cosbval = cos(b * degtorad)
            if cosbval * coslval < costhetaval:
                incl = 1. - (1. - pow(rcut/rodot,2) ) * pow(coslval*cosbval,2)
                if incl > 0:
                    integral = 0.0
                    s = smin
                    for l in range(Ns):

                        flux_min_eval = max(Lmin/s**2,4*pi*fluxmin)
                        flux_max_eval = min(Lmax/s**2,4*pi*fluxmax)
                        if flux_max_eval < flux_min_eval:
                            pref_L = 0.0
                        else:
                            pref_L = (pow(flux_min_eval, 1.-beta) - pow(flux_max_eval, 1.-beta))*pref_L_norm

                        r_squared = (s*cosbval*coslval - rodot)**2 + s**2*(1. - pow(cosbval*coslval,2) )
                        if r_squared < rcut**2:
                            integral += pref_L * pow(s,4.-2.*beta)*pow(r_squared,-alpha/2.)
                        s += ds

                    total_res += cosbval * integral 

    return total_res * d_ang**2. * degtorad**2 * ds * pref_rho







######################################################### 
# Predicted Disk PSRs for prior                         #
#########################################################
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef inline double Ndisk_total(double Ndisk, double n, double sigma,double z0, double beta, double Lmin, double Lmax ,int Ns ,int Nang ,double smax, double theta_mask ) nogil:
    """ Return the number of bulge PSRs in bin k for given parameters
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
    cdef double fluxmin = 1.8e-5 * fluxunits #MeV cm^{-2} s^{-1}
    cdef double fluxmax = 0.003896 * fluxunits  #MeV cm^{-2} s^{-1} #6.539e-4

    # Setup output and if variables
    cdef double Nijk = 0.
    cdef double prefactor, pref_rho, pref_L, lim_m, lim_p, integral,s,ds

    cdef int l

    cdef double res = 0.0

    cdef int i_b, i_ell
    cdef double ell,b
    cdef double ell_start = -180.0 #integrate from -180 to +180 in lon
    cdef double b_start = -90.0 #and from -90 to +90 in lat
    cdef double d_ang_ell = (360.0)/float(Nang)
    cdef double d_ang_b = (180.0)/float(Nang)

    cdef double total_res = 0.0

    cdef double coslval, cosbval, sinbval

    ####
    ##cos theta for the mask
    cdef double costhetaval = cos(theta_mask * degtorad)

    cdef double R,z


    # Common prefactors
    pref_rho = Ndisk /4./pi/z0/pow(sigma,n+2)/tgamma(n+2)
    #pref_L = pow(4.*pi,1.-beta)*(pow(fluxmin, 1.-beta) - pow(fluxmax, 1.-beta))/(pow(Lmin,1-beta) - pow(Lmax,1-beta)) #(beta-1.)
    cdef double pref_L_norm = 1./(pow(Lmin,1-beta) - pow(Lmax,1-beta))

    cdef double flux_min_eval, flux_max_eval

    # for s integration
    smin = 0.01
    ds = (smax - smin)/float(Ns)

    for i_b in range(Nang):
        b = b_start + i_b*d_ang_b + d_ang_b/2.
        for i_ell in range(Nang):
            ell = ell_start + i_ell*d_ang_ell + d_ang_ell/2.
            coslval = cos(ell * degtorad)
            cosbval = cos(b * degtorad)
            if cosbval * coslval < costhetaval:

                res = 0.0
                integral = 0.0
                s = smin
                for l in range(Ns):

                    flux_min_eval = max(Lmin/s**2,4*pi*fluxmin)
                    flux_max_eval = min(Lmax/s**2,4*pi*fluxmax)
                    if flux_max_eval < flux_min_eval:
                        pref_L = 0.0
                    else:
                        pref_L = (pow(flux_min_eval, 1.-beta) - pow(flux_max_eval, 1.-beta))*pref_L_norm

                    r_squared = (s*cosbval*coslval - rodot)**2 + s**2*(1. - pow(cosbval*coslval,2) )
                    z = s*sqrt(1 - cosbval**2)
                    R = sqrt( r_squared - z**2 )
                    integral += pref_L * pow(s,4.-2.*beta)*pow(R,n)*exp(-R/sigma - z/z0)
                    s += ds

                res = cosbval * integral 
                total_res += cosbval * integral

    return total_res * d_ang_ell*d_ang_b * degtorad**2 * ds  * pref_rho


