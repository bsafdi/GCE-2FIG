###############################################################################
# bin_data.py
###############################################################################
#
# Bin the 2FIG PSR candidates into the l, b and flux bins
#
# Written: Nick Rodd, MIT, 5 August 2017
#
###############################################################################

import numpy as np

# Set mask in degrees - distance from GC within which we ignore events
mask = 2.

# Load the raw data, contains 86 pulsar candidates with associated l, b, flux values
rd = np.loadtxt('./psrcandidates.dat')

# Define our bin edges
angbins = -np.linspace(-20,20,13) # - because top left is 0,0 and astro longitude
fluxbins = np.append(np.logspace(-6,-5,7), np.logspace(-5,-4,3)[1:3])

PSR_data = np.zeros((12,12,8))

# Now bin the data
for ips in range(len(rd)):
    lval = rd[ips, 0]
    bval = rd[ips, 1]
    fval = rd[ips, 2]
    if np.cos(lval*np.pi/180.)*np.cos(bval*np.pi/180.) > np.cos(mask*np.pi/180.): continue
    for li in range(12):
        if (lval <= angbins[li]) & (lval > angbins[li+1]):
            for bi in range(12):
                if (bval <= angbins[bi]) & (bval > angbins[bi+1]):
                    for fi in range(8):
                        if (fval >= fluxbins[fi]) & (fval < fluxbins[fi+1]):
                            PSR_data[li,bi,fi] += 1.

np.save('./PSR_data',PSR_data)
