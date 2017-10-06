###############################################################################
# PSR_flux.pyx
###############################################################################
#
# Return the predicted bulge flux in bin (i,j,k).
# Analogous to PSR_counts.pyx
#
# Data is binned from top left (l,b)=(+20,+20) to bottom right (-20,-20)
#
###############################################################################


# Import basic functions
import numpy as np
cimport numpy as np
cimport cython

# C functions
cdef extern from "math.h":
    double pow(double x, double y) nogil
    double sqrt(double x) nogil
    double cos(double x) nogil
    double sin(double x) nogil
    double tgamma(double x) nogil
    double exp(double x) nogil
    double fabs(double x) nogil

# Useful variables
cdef double pi = np.pi
cdef double rodot = 8.5 # Earth GC distance in kpc
cdef double degtorad = pi/180.
cdef double fluxunits = 1.52529e37 # [converts MeV/cm^2 to erg/kpc^2]

# Angular bin boundaries for b and l in [degrees]
# +20 to -20 because we work from the top left to the bottom right
cdef double[::1] ang_boundaries = -np.linspace(-20.,20.,13,dtype=np.float)

# Flux bin boundaries in [Mev/cm^2/s]
cdef double fluxvals[9]
fluxvals = [1.00000000e-06, 1.46779927e-06, 2.15443469e-06, 3.16227766e-06, 
            4.64158883e-06, 6.81292069e-06, 1.00000000e-05, 3.16227766e-05, 
            1.00000000e-04]


#########################################################
# Predicted Bulge PSRs in bin (i,j,k) [long, lat, flux] #
#########################################################

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef inline double Fbulge_ijk(int i, int j, int k, double Nbulge,
                               double omega_ijk, double alpha, double beta,
                               double rcut, double Lmin, double Lmax, int Ns,
                               int Nang, double theta_mask) nogil:
    """ Return the number of bulge PSRs in bin (i,j,k) for given parameters
    
    ---arguments---
    Ns : number of steps in s integration
    Nang : number of steps in the angular integration
    
    ----units---
    rcut : in kpc, distance from GC after which to set bulge to 0
    Lmin : [erg/s]
    Lmax : [erg/s]
    theta_mask : degrees, mask from the GC   
    """
    
    ## Setup variables
    cdef double fluxmin = fluxvals[k] * fluxunits # min flux in [erg/kpc^2/s]
    cdef double fluxmax = fluxvals[k+1] * fluxunits # max flux in [erg/kpc^2/s]

    cdef double l_start = ang_boundaries[i] # initial l value
    cdef double b_start = ang_boundaries[j] # initial b value
    cdef double dang = (ang_boundaries[0] - ang_boundaries[1])/float(Nang)
    cdef double l, b, coslval, cosbval

    cdef double smin = rodot - rcut # min possible s value for LOS integral
    cdef double smax = rodot + rcut # max possible s value for LOS integral
    cdef double ds = (smax - smin)/float(Ns)

    cdef double cosmask = cos(theta_mask * degtorad) # mask cos theta

    cdef Py_ssize_t bi, li, si # loop indexes
    cdef double inbulge, sum_1ang, fmin_s, fmax_s, pref_f, rsq

    cdef double sum_tot = 0. # running total of the sum to get the integral

    # Prefactors for the spatial and luminosity distributions
    # RTB: changed power of beta -> beta-1
    cdef double pref_rho = omega_ijk * Nbulge * (3.-alpha) \
                           / pow(4.*pi, beta-1.) / pow(rcut, 3.-alpha)
    
    # RTB: defined to integrate to unity for number of sources! (1-beta drops out for number of sources)
    cdef double pref_L = 1./(pow(Lmax, 1.-beta) - pow(Lmin, 1.-beta))


    ## Perform the integral as a discretized sum in l, b, and s
    for li in range(Nang):
        l = l_start - li*dang - dang/2.
        coslval = cos(l * degtorad)
        for bi in range(Nang):
            b = b_start - bi*dang - dang/2.
            cosbval = cos(b * degtorad)

            # Only proceed if outside the mask (so a smaller cosine)
            if cosbval * coslval > cosmask: continue

            # Only proceed if this angle has some overlap with the bulge
            inbulge = pow(cosbval*coslval, 2.) + pow(rcut/rodot, 2.) - 1.
            if inbulge <= 0.: continue
                
            # Perform the sum over s
            sum_1ang = 0.
            sval = smin
            for si in range(Ns):
                # Determine appropriate flux boundaries at this distance,
                # ensuring we remain within our flux boundaries
                fmin_s = max(Lmin/pow(sval,2.)/4./pi, fluxmin)
                fmax_s = min(Lmax/pow(sval,2.)/4./pi, fluxmax)

                # Check if this s value puts us outside a sensible flux range
                if fmax_s < fmin_s:
                    pref_f = 0.0
                else:
                    # RTB: changed powers and added additional factors
                    pref_f = (pow(fmax_s, 2.-beta) - pow(fmin_s, 2.-beta)) * pref_L * (1. - beta) / (2. - beta) / (4. * pi * pow(sval,2.))

                # Determine distance squared from GC for this point
                rsq = pow(sval,2.) - 2.*sval*rodot*cosbval*coslval + pow(rodot,2.)
                
                # This point only contributes if it is within the bulge
                if rsq < pow(rcut, 2.):
                    # RTB: added in extra factor of s^2 here
                    sum_1ang += pref_f * pow(sval,6.-2.*beta)*pow(rsq,-alpha/2.)
                
                sval += ds # increment s for the next step

            sum_tot += cosbval * sum_1ang # account for cosb in Jacobian 

    # Return the integral - sum times the step sizes and the prefactor
    # RTB: convert units back again.
    return sum_tot * pow(dang,2.) * pow(degtorad,2.) * ds * pref_rho / fluxunits