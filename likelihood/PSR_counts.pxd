###############################################################################
# PSR_counts.pxd
###############################################################################
#
# Here we predefine all functions in PSR_counts.pyx so that they are compiled
# simultaneously - this way the code optimizes all functions at once and allows
# functions to be called as pure C
#
###############################################################################

cpdef double Nbulge_ijk(int, int, int, double, double, double, double, double, 
                        double, double,int,int, double) nogil
    
cpdef double Ndisk_ijk(int, int, int, double, double, double, double, double, 
                       double, double, double, int, int, double, double) nogil

cpdef double Nbulge_total(double, double, double, double, double, double, int, 
                          int, double) nogil

cpdef double Ndisk_total(double, double, double, double, double, double, double, 
                         int, int, double, double) nogil
