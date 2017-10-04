import sys,os
import numpy as np
import matplotlib
matplotlib.use('Agg')

from run import run_scan

fixed_params = ['n','sigma','Lmax_disk', 'Lmax_bulge']
fixed_param_vals = [2.35,1.528,1.0e36,1.0e36]

floated_params = ['N_bulge','N_disk', 'alpha','z0', 'beta']
floated_param_priors = [[0,3000000],[0,3000000],[2.1,5.0],[0.01,2.0],[1.1,3.0]]

rs_nd_nb_alpha = run_scan(fixed_params, fixed_param_vals, floated_params, floated_param_priors, Ns = 500, Nang = 10, share_betas=True, use_prior=True, Nang_prior=500)
rs_nd_nb_alpha.perform_scan_multinest(nlive=1500, chains_dir='chains/rs_nd_nb_alpha_nlive1500/')
