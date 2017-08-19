import sys,os
import numpy as np
import matplotlib
matplotlib.use('Agg')

from run import run_scan

fixed_params = ['n','sigma','alpha','Lmax_disk', 'Lmax_bulge','N_bulge']
fixed_param_vals = [2.35,1.528,2.6,1.0e36,1.0e36,0]

floated_params = ['N_disk', 'z0', 'beta']
floated_param_priors = [[0,2000000],[0.01,2.0],[1.1,3.0]]

rs_nd = run_scan(fixed_params, fixed_param_vals, floated_params, floated_param_priors, Ns = 100, Nang = 4, share_betas=True, use_prior=True)

rs_nd.perform_scan_minuit()

fixed_params = ['n','sigma','alpha','Lmax_disk', 'Lmax_bulge']
fixed_param_vals = [2.35,1.528,2.6,1.0e36,1.0e36]

floated_params = ['N_bulge','N_disk','z0', 'beta']
floated_param_priors = [[0,2000000],[0,2000000],[0.01,2.0],[1.1,3.0]]

rs_nd_nb = run_scan(fixed_params, fixed_param_vals, floated_params, floated_param_priors, Ns = 100, Nang = 4, share_betas=True, use_prior=True)

rs_nd_nb.perform_scan_minuit()

TS = 2*(rs_nd_nb.max_LL - rs_nd.max_LL)
print "TS = ", TS
print np.sqrt(TS), "sigma"

fixed_params = ['n','sigma','Lmax_disk', 'Lmax_bulge']
fixed_param_vals = [2.35,1.528,1.0e36,1.0e36]

floated_params = ['N_bulge','N_disk', 'alpha','z0', 'beta']
floated_param_priors = [[0,3000000],[0,3000000],[2.1,3.0],[0.01,2.0],[1.1,2.0]]

rs_nd_nb_alpha = run_scan(fixed_params, fixed_param_vals, floated_params, floated_param_priors, Ns = 100, Nang = 4, share_betas=True, use_prior=True)

rs_nd_nb_alpha.perform_scan_minuit()

TS = 2*(rs_nd_nb_alpha.max_LL - rs_nd.max_LL)
print "TS = ", TS
print np.sqrt(TS), "sigma"