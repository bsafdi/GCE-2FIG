###############################################################################
# PSR_base_analysis.py
###############################################################################
#
# Perform the analysis using the likelihood coded up following 1705.00009
# Here we only scan over Nbulge and Ndisk, leaving other paramters fixed
#
# Written: Nick Rodd, MIT, 5 August 2017
#
###############################################################################

# Fix base parameters
n=2.35
sigma=1.528
alpha=2.6
z0=0.7
beta=1.2

# Set up log priors on Ndisk and Nbulge
Nbulge_prior = [0,1e7]
Ndisk_prior = [0,1e7]


# Basic modules
import numpy as np
import pymultinest

# Import custom ll
import sys
sys.path.append('/group/hepheno/smsharma/GCE-2FIG-Ben/python')
import LL

# Load the binned efficiency and data files
PSR_data = np.load('../data/PSR_data.npy') + np.load('../data/PSR_data_3fgl.npy')
omega_ijk = np.load('../data/omega_ijk.npy')

# Multinest settings
nlive = 200
chains_dir = './chains/base/'
pymultinest_options = {'importance_nested_sampling': False,
                        'resume': False, 'verbose': True,
                        'sampling_efficiency': 'model',
                        'init_MPI': False, 'evidence_tolerance': 0.5,
                        'const_efficiency_mode': False}

# Setup priors
theta_min = [Nbulge_prior[0], Ndisk_prior[0]]
theta_max = [Nbulge_prior[1], Ndisk_prior[1]]
theta_interval = list(np.array(theta_max) - np.array(theta_min))

def prior_cube(cube, ndim=1, nparams=1):
    """ Cube of priors - in the format required by MultiNest
    """

    for i in range(ndim):
        cube[i] = cube[i] * theta_interval[i] + theta_min[i]
    return cube

# Define likelihood
def ll(theta, ndim = 1, nparams = 1):
    """ Log Likelihood in the format required by Multinest
    """

    return LL.ll(PSR_data, omega_ijk, theta[0], theta[1],
                     alpha, n, sigma, z0, beta, beta)

# Run multinest
n_params = len(theta_min)
pymultinest.run(ll, prior_cube, n_params, outputfiles_basename = chains_dir, 
                n_live_points = nlive, **pymultinest_options)