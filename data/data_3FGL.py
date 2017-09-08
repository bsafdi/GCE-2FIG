###############################################################################
# data_3FGL.py
###############################################################################
#
# Bin 3FGL psr candidates using file from Mattia
#
# Data is binned from top left (l,b)=(+20,+20) to bottom right (-20,-20)
#
# Written: Nick Rodd, MIT, 7 August 2017
#
###############################################################################

import numpy as np
from astropy.io import fits

# Set mask in degrees - distance from GC within which we ignore events
mask = 2.

# Load Mattia catalog
load = fits.open('./2FIG_Pass8_Arxiv.fits')
cat = np.array(load[1].data)

glon = np.array([])
glat = np.array([])
eflux = np.array([])

count = 0

for i in range(len(cat)):
    # Check if PSR or psr
    if (cat[i][22] == 'PSR') | (cat[i][22] == 'psr'):
        glon = np.append(glon,cat[i][3])
        glat = np.append(glat,cat[i][4])
        eflux = np.append(eflux,cat[i][11])
        count += 1

# Now bin
angbins = -np.linspace(-20,20,13) # - because top left is 0,0 and astro longitude
fluxbins = np.append(np.logspace(-6,-5,7), np.logspace(-5,-4,3)[1:3])

PSR_data = np.zeros((12,12,8))

for ips in range(len(glon)):
    lval = glon[ips]
    bval = glat[ips]
    if np.cos(lval*np.pi/180.)*np.cos(bval*np.pi/180.) > np.cos(mask*np.pi/180.): continue
    fval = eflux[ips]
    for li in range(12):
        if (lval <= angbins[li]) & (lval > angbins[li+1]):
            for bi in range(12):
                if (bval <= angbins[bi]) & (bval > angbins[bi+1]):
                    for fi in range(8):
                        if (fval >= fluxbins[fi]) & (fval < fluxbins[fi+1]):
                            PSR_data[li,bi,fi] += 1.

np.save('./PSR_data_3fgl',PSR_data)
