###############################################################################
# get_counts_inline.pxd
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

cpdef double Nbulge_full_ang_ijk(int , int , int , double , double , double , double ,double  , double , double  ,int  ,int , double ) nogil
    
cpdef double Ndisk_full_ang_ijk(int , int , int , double , double , double, double ,double , double , double , double ,int ,int ,double , double  ) nogil

cpdef double Nbulge_total(double , double , double ,double , double , double ,int ,int , double) nogil

cpdef double Ndisk_total(double, double, double,double, double, double, double ,int,int ,double, double) nogil
