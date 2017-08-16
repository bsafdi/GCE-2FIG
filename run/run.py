###############################################################################
# run.py
###############################################################################
#
# Perform scan to obtain best-fit disk and bulge PS population parameters
#
# Written: Siddharth Mishra-Sharma, Princeton, 15 August 2017
# 
###############################################################################

import sys, os
sys.path.append("../python/")

import numpy as np
import pymultinest

import LL

class run_scan():
    def __init__(self, fixed_params, fixed_param_vals, floated_params, floated_param_priors, data_dir = '../data/', chains_dir = '/group/hepheno/smsharma/GCE-2FIG-bsafdi/run/chains/analysis/',
                 Lmin = 1.0e31, Lmax_disk = 1.0e36, Lmax_bulge = 1.0e36, Ns = 100, Nang = 1, smax_disk = 40, theta_mask=2.0):
        """ Initialize scan class.

            :param fixed_params: Array of parameters to be held fixed
            :param fixed_param_vals: Array of values for parameters to be held fixed
            :param floated_params: Array of parameters to be floated
            :param floated_param_priors: Priors array for parameters to be floated
        """
        self.fixed_params = fixed_params
        self.fixed_param_vals = fixed_param_vals
        self.floated_params = floated_params
        self.floated_param_priors = floated_param_priors
        self.data_dir = data_dir

        self.Lmin = Lmin
        self.Lmax_disk = Lmax_disk
        self.Lmax_bulge = Lmax_bulge
        self.Ns = Ns
        self.Nang = Nang
        self.smax_disk = smax_disk
        self.theta_mask = theta_mask
        self.chains_dir = chains_dir

        self.all_params = ['N_bulge','N_disk','alpha','n','sigma','z0','beta_disk','beta_bulge']

        self.make_dirs([self.chains_dir])
        self.load_data()
        self.setup_fixed_params()
        self.setup_prior_thetas()

    def load_data(self):
        """ Load the binned efficiency and data files
        """       
        self.PSR_data = np.load(self.data_dir + '/PSR_data.npy') + np.load(self.data_dir + '/PSR_data_3fgl.npy')
        self.omega_ijk = np.load(self.data_dir + '/omega_ijk.npy')

    def setup_fixed_params(self):
        """ Set up values and arrays for parameters to be held fixed
        """ 
        self.theta = np.zeros(len(self.all_params))
        for param in self.all_params:
            if param in self.fixed_params:
                param_pos = np.argwhere(np.array(self.fixed_params) == param)[0][0]
                param_all_pos = np.argwhere(np.array(self.all_params) == param)[0][0]
                self.theta[param_all_pos] = self.fixed_param_vals[param_pos]

    def setup_prior_thetas(self):
        """ Set up priors and arrays for parameters to be floated
        """ 
        self.param_all_pos_ary = []
        self.theta_min = []
        self.theta_max = []
        for param in self.all_params:
            if param in self.floated_params:
                param_pos = np.argwhere(np.array(self.floated_params) == param)[0][0]
                param_all_pos = np.argwhere(np.array(self.all_params) == param)[0][0]
                self.param_all_pos_ary.append(param_all_pos)
                self.theta_min += [self.floated_param_priors[param_pos][0]]
                self.theta_max += [self.floated_param_priors[param_pos][1]]

        self.theta_interval = list(np.array(self.theta_max) - np.array(self.theta_min))

    def prior_cube(self, cube, ndim=1, nparams=1):
        """ Cube of priors - in the format required by MultiNest
        """
        for i in range(ndim):
            cube[i] = cube[i] * self.theta_interval[i] + self.theta_min[i]
        return cube

    def ll(self, theta, ndim = 1, nparams = 1):
        """ Log Likelihood in the format required by MultiNest
        """
        theta_ll = self.theta
        for float_idx, float_item in enumerate(self.param_all_pos_ary):
            theta_ll[float_item] = theta[float_idx]
        
        ll_val =  LL.ll(self.PSR_data, self.omega_ijk,                  
                        *theta_ll, 
                        Lmin = self.Lmin, Lmax_disk = self.Lmax_disk, Lmax_bulge = self.Lmax_bulge, 
                        Ns = self.Ns, Nang = self.Nang, smax_disk = self.smax_disk, theta_mask = self.theta_mask)
        print theta_ll, ll_val

        return ll_val

    def perform_scan(self):
        # Run multinest
        n_params = len(self.floated_params)
        nlive = 200
        pymultinest_options = {'importance_nested_sampling': False,
                                'resume': False, 'verbose': True,
                                'sampling_efficiency': 'model',
                                'init_MPI': False, 'evidence_tolerance': 0.5,
                                'const_efficiency_mode': False}

        pymultinest.run(self.ll, self.prior_cube, n_params, outputfiles_basename = self.chains_dir, 
                        n_live_points = nlive, **pymultinest_options)

    @staticmethod
    def make_dirs(dirs):
        """ Creates directories if they do not already exist 
        """
        for d in dirs:
            if not os.path.exists(d):
                try:
                    os.mkdir(d)
                except OSError as e:
                    if e.errno != 17:
                        raise

