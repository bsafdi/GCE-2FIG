import sys,os
import numpy as np
import matplotlib
matplotlib.use('Agg')

from run import run_scan

fixed_params = ['n','sigma','alpha','Lmax_disk', 'Lmax_bulge','N_bulge']
fixed_param_vals = [2.35,1.528,2.6,1.0e36,1.0e36,0]

floated_params = ['N_disk', 'z0', 'beta']
floated_param_priors = [[0,3000000],[0.01,2.0],[1.1,3.0]]

rs_nd = run_scan(fixed_params, fixed_param_vals, floated_params, floated_param_priors, Ns = 500, Nang = 10, share_betas=True, use_prior=True, Nang_prior=500, efficiency_file = "efficiency_long_dan_trunc.npy")
rs_nd.perform_scan_multinest(nlive=1200, chains_dir='chains/rs_nd_2PCDan/')
