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

#########################################################
# This is New Code written by Ben Safdi                 # 
# Predicted Bulge PSRs in bin (i,j,k) [long, lat, flux] #
# Performs the angular integral instead of taking the   #
# value from the bin centre                             #
#########################################################
cpdef double Nbulge_full_ang_ijk(int , int , int , double , double , double , double ,double  , double , double  ,int  ,int , double ) nogil
    



#########################################################
# This is New Code written by Ben Safdi                 # 
# Predicted Disk PSRs in bin (i,j,k) [long, lat, flux] #
# Performs the angular integral instead of taking the   #
# value from the bin centre                             #
#########################################################
cpdef double Ndisk_full_ang_ijk(int , int , int , double , double , double, double ,double , double , double , double ,int ,int ,double , double  ) nogil


cpdef double Nbulge_total(double , double , double ,double , double , double ,int ,int , double) nogil


cpdef double Ndisk_total(double, double, double,double, double, double, double ,int,int ,double, double) nogil








