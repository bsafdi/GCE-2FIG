###############################################################################
# bulge_pscounts.pyx
###############################################################################
#
# Return the predicted bulge counts in a bin (i,j,k), following the definition
# in 1705.00009
#
# Written: Nick Rodd, MIT, 5 August 2017
# Moddified by Ben Safdi, UM August 15 2017
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
#########################################################
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef double Nbulge_full_ijk(int i, int j, int k, double Nbulge, double omega_ijk, double alpha, double beta,double rcut = 3.0, double Lmin= 1.0e31, double Lmax = 1.0e36,int Ns = 1000) nogil:
    """ Return the number of bulge PSRs in bin (i,j,k) for given parameters
    ---arguments---
    Ns : number of points in s integration
    ----units---
    rcut : in kpc
    Lmin : erg s^{-1}
    Lmax : erg s^{-1}
    """
    # Determine angles and flux boundaries, convert flux to erg/kpc^2/s
    cdef double coslval = cos(angvals[i] * degtorad)
    cdef double cosbval = cos(angvals[j] * degtorad)
    cdef double fluxmin = fluxvals[k] * fluxunits
    cdef double fluxmax = fluxvals[k+1] * fluxunits

    # Setup output and if variables
    cdef double Nijk = 0.
    cdef double prefactor, pref_rho, pref_L, lim_m, lim_p, integral, smin, smax,s,ds

    cdef int l

    cdef double incl = 1. - (1. - pow(rcut/rodot,2) ) * pow(coslval*cosbval,2)

    if incl > 0:
        pref_rho = omega_ijk * Nbulge * (3. - alpha) / 4. / pi / pow(rcut,3 - alpha)
        pref_L = pow(4.*pi,1.-beta)*(pow(fluxmin, 1.-beta) - pow(fluxmax, 1.-beta))/(pow(Lmin,1-beta) - pow(Lmax,1-beta)) #(beta-1.)

        smin = rodot - rcut
        smax = rodot + rcut
        ds = (smax - smin)/float(Ns)
        integral = 0.0
        s = smin
        for l in range(Ns):
            r_squared = (s*cosbval*coslval - rodot)**2 + s**2*(1. - pow(cosbval*coslval,2) )
            #print r_squared
            if r_squared < rcut**2:
                integral += pow(s,4.-2.*beta)*pow(r_squared,-alpha/2.)
            s += ds

        res = area * cosbval * integral * ds * pref_L * pref_rho
        #print integral, ds, pref_L, pref_rho

    #print area, cosbval, integral, ds, pref_L, pref_rho
    return res


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
cpdef double Nbulge_full_ang_ijk(int i, int j, int k, double Nbulge, double omega_ijk, double alpha, double beta,double rcut = 3.0, double Lmin= 1.0e31, double Lmax = 1.0e36,int Ns = 1000,int Nang = 10) nogil:
    """ Return the number of bulge PSRs in bin (i,j,k) for given parameters
    ---arguments---
    Ns : number of points in s integration
    ----units---
    rcut : in kpc
    Lmin : erg s^{-1}
    Lmax : erg s^{-1}
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


    # Common prefactors
    pref_rho = omega_ijk * Nbulge * (3. - alpha) / 4. / pi / pow(rcut,3 - alpha)
    pref_L = pow(4.*pi,1.-beta)*(pow(fluxmin, 1.-beta) - pow(fluxmax, 1.-beta))/(pow(Lmin,1-beta) - pow(Lmax,1-beta)) #(beta-1.)

    # for s integration
    smin = rodot - rcut
    smax = rodot + rcut
    ds = (smax - smin)/float(Ns)

    for i_b in range(Nang):
        b = b_start + i_b*d_ang
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
cpdef double Ndisk_full_ang_ijk(int i, int j, int k, double Ndisk, double omega_ijk, double n, double sigma,double z0, double beta, double Lmin= 1.0e31, double Lmax = 1.0e36,int Ns = 1000,int Nang = 10,double smax=40) nogil:
    """ Return the number of bulge PSRs in bin (i,j,k) for given parameters
    ---arguments---
    Ns : number of points in s integration
    ----units---
    z0 : kpc
    sigma : kpc
    Lmin : erg s^{-1}
    Lmax : erg s^{-1}
    smax : 40 kpc default, how far to integrate out to from Sun?
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

    cdef double R,z


    # Common prefactors
    pref_rho = omega_ijk * Ndisk /4./pi/z0/pow(sigma,n+2)/tgamma(n+2)
    pref_L = pow(4.*pi,1.-beta)*(pow(fluxmin, 1.-beta) - pow(fluxmax, 1.-beta))/(pow(Lmin,1-beta) - pow(Lmax,1-beta)) #(beta-1.)

    # for s integration
    smin = 0.0
    ds = (smax - smin)/float(Ns)

    for i_b in range(Nang):
        b = b_start + i_b*d_ang
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

























#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
# Old code from Nick Below                              #


#########################################################
# Predicted Bulge PSRs in bin (i,j,k) [long, lat, flux] #
#########################################################

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef double Nbulge_ijk(int i, int j, int k, double Nbulge, double omega_ijk,
                        double alpha, double beta) nogil:
    """ Return the number of bulge PSRs in bin (i,j,k) for given parameters
    """
    
    # Determine angles and flux boundaries, convert flux to erg/kpc^2/s
    cdef double coslval = cos(angvals[i] * degtorad)
    cdef double cosbval = cos(angvals[j] * degtorad)
    cdef double fluxmin = fluxvals[k] * fluxunits
    cdef double fluxmax = fluxvals[k+1] * fluxunits

    # Setup output and if variables
    cdef double Nijk = 0.
    cdef double prefactor, lim_m, lim_p, integral

    # Figure out if line of sight includes bulge, if not return 0
    cdef double sdet = pow(coslval * cosbval, 2.) + 9. / pow(rodot, 2.) \
                       - 1.
    if sdet > 0.:        
        # Compute prefactor
        prefactor = Nbulge * omega_ijk * 100. * (3.-alpha) \
                    / pow(4*pi, beta) / pow(3., 5.-alpha) \
                    * (pow(fluxmax, 1.-beta) - pow(fluxmin, 1.-beta)) \
                    / (pow(1.0e36, 1.-beta) - pow(1.0e31, 1.-beta))

        # Compute limits and integral
        lim_m = rodot * (coslval * cosbval - sqrt(sdet))
        lim_p = rodot * (coslval * cosbval + sqrt(sdet))

        integral = bulge_int(alpha, beta, coslval, cosbval, lim_m, lim_p)

        # Combine factors and convert units to sr from degrees^2
        Nijk = prefactor * integral * pow(degtorad,2.)

    return Nijk


########################################################
# Predicted Disk PSRs in bin (i,j,k) [long, lat, flux] #
########################################################

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef double Ndisk_ijk(int i, int j, int k, double Ndisk, double omega_ijk,
                       double n, double sigma, double z0, double beta) nogil:
    """ Return the number of disk PSRs in bin (i,j,k) for given parameters
    """
    
    # Determine angles and flux boundaries, convert flux to erg/kpc^2/s
    cdef double coslval = cos(angvals[i] * degtorad)
    cdef double cosbval = cos(angvals[j] * degtorad)
    cdef double tanbval = fabs(tan(angvals[j] * degtorad)) # abs as in |z|
    cdef double fluxmin = fluxvals[k] * fluxunits
    cdef double fluxmax = fluxvals[k+1] * fluxunits

    # Compute prefactor
    cdef double prefactor = Ndisk * omega_ijk * 100. / (9. * pow(4*pi, beta)
                            * z0 * pow(sigma, 2.+n) * tgamma(2.+n)) \
                            * (pow(fluxmax, 1.-beta) - pow(fluxmin, 1.-beta)) \
                            / (pow(1.0e36, 1.-beta) - pow(1.0e31, 1.-beta))

    # Compute integral
    cdef double integral = disk_int(n, sigma, z0, beta, coslval, cosbval, tanbval) 

    # Combine factors and convert units to sr from degrees^2
    cdef double Nijk = prefactor * integral * pow(degtorad,2.)

    return Nijk


###########################
# Calculation of Integral #
###########################

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double bulge_int(double alpha, double beta, double coslval, double cosbval,
                      double lim_m, double lim_p) nogil:
    """ Calculate the integral over the bulge from lim_m to lim_p
        For speed do the integral by the rectangle rule
    """
    
    cdef double intval = 0.
    cdef double sval = lim_m + sbin/2.
    while (sval < lim_p):
        intval += integrand(sval, alpha, beta, coslval, cosbval)
        sval += sbin

    intval *= sbin

    return intval

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double integrand(double sval, double alpha, double beta, double coslval, 
                      double cosbval) nogil:
    """ Evaluate the bulge integrand
    """

    return pow(sval, 4.-2.*beta) \
           * pow(pow(sval, 2.) - 2*rodot*sval*coslval*cosbval + pow(rodot, 2.),
                 -alpha/2.)



###########################
# Calculation of Integral #
###########################

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double disk_int(double n, double sigma, double z0, double beta, 
                     double coslval, double cosbval, double tanbval) nogil:
    """ Calculate the integral over the disk from 0 to smax
        For speed do the integral by the rectangle rule
    """
    
    cdef double intval = 0.
    cdef double sval = sbin/2.
    cdef double smax = 30.0
    while (sval < smax):
        intval += integrand_disk(sval, n, sigma, z0, beta, coslval, cosbval, tanbval)
        sval += sbin

    intval *= sbin

    return intval

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double integrand_disk(double sval, double n, double sigma, double z0, 
                      double beta, double coslval, double cosbval, 
                      double tanbval) nogil:
    """ Evaluate the disk integrand
    """

    cdef double Rval = sqrt(pow(sval, 2.) - 2*rodot*sval*coslval + pow(rodot, 2.))
    cdef double zval = sval * tanbval

    return pow(sval, 4.-2.*beta) * pow(Rval, n) * exp(-Rval/sigma) * exp(-zval/z0) 