# Extract data from Sid runs

import sys,os
import numpy as np
import matplotlib
matplotlib.use('Agg')

from run import run_scan

# Load the no bulge case
fixed_params = ['n','sigma','alpha','Lmax_disk', 'Lmax_bulge','N_bulge']
fixed_param_vals = [2.35,1.528,2.6,1.e36,1.e36,0.]

floated_params = ['N_disk', 'z0', 'beta']
floated_param_priors = [[0.,1.e7],[0.01,2.],[1.1,3.]]

rs_nd = run_scan(fixed_params, fixed_param_vals, floated_params, floated_param_priors, Ns = 200, Nang = 5, share_betas=True, use_prior=True)

print rs_nd.get_max_log_likelihood(chains_dir='chains/np_nb_nlive500/')

LL_nb = rs_nd.get_max_log_likelihood(chains_dir='chains/np_nb_nlive500/')['log_likelihood']

print "LL_nb:",2*LL_nb

# Load the bulge case
fixed_params = ['n','sigma','alpha','Lmax_disk', 'Lmax_bulge']
fixed_param_vals = [2.35,1.528,2.6,1.e36,1.e36]

floated_params = ['N_bulge','N_disk', 'z0', 'beta']
floated_param_priors = [[0.,1.e7],[0.,1.e7],[0.01,2.],[1.1,3.]]

rs_wb = run_scan(fixed_params, fixed_param_vals, floated_params, floated_param_priors, Ns = 200, Nang = 5, share_betas=True, use_prior=False)

print rs_wb.get_max_log_likelihood(chains_dir='chains/np_nlive500/')

LL_wb = rs_wb.get_max_log_likelihood(chains_dir='chains/np_nlive500/')['log_likelihood']

print "LL_wb:",2*LL_wb

print "TS:",2*(LL_wb-LL_nb)
