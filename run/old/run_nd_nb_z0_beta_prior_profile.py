import sys,os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')

from run import run_scan

parser = argparse.ArgumentParser()
parser.add_argument("--Nb",
                  action="store", dest="Nb", default=1.0,type=float)

results = parser.parse_args()
Nb=results.Nb

fixed_params = ['N_bulge','n','sigma','alpha','Lmax_disk', 'Lmax_bulge']
fixed_param_vals = [Nb, 2.35,1.528,2.6,1.0e36,1.0e36]

floated_params = ['N_disk','z0', 'beta']
floated_param_priors = [[0,3000000],[0.01,2.0],[1.1,3.0]]

rs_nd_nb = run_scan(fixed_params, fixed_param_vals, floated_params, floated_param_priors, Ns = 1000, Nang = 10, share_betas=True, use_prior=True, Nang_prior=500)
rs_nd_nb.perform_scan_multinest(nlive=500, chains_dir='chains/rs_nd_nb_morepoints_nb_' + str(int(Nb))[:10]+'/')
