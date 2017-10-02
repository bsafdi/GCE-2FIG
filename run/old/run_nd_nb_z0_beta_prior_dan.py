import sys,os
import numpy as np
import matplotlib
matplotlib.use('Agg')

from run import run_scan

fixed_params = ['n','sigma','alpha','Lmax_disk', 'Lmax_bulge']
fixed_param_vals = [2.35,1.528,2.6,1.0e36,1.0e36]

floated_params = ['N_bulge','N_disk','z0', 'beta']
floated_param_priors = [[0,3000000],[0,3000000],[0.01,2.0],[1.1,3.0]]

rs_nd_nb = run_scan(fixed_params, fixed_param_vals, floated_params, floated_param_priors, Ns = 1000, Nang = 10, share_betas=True, use_prior=True, Nang_prior=500, efficiency_file = "efficiency_long_dan.npy")
rs_nd_nb.perform_scan_multinest(nlive=500, chains_dir='chains/rs_nd_nb_morepoints_2PCeff_Dan/')
