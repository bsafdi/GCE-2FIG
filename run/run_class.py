# Basic modules
import numpy as np
import pymultinest
import sys
sys.path.append("../python/")
import get_counts_inline as gc
import LL

class run_scan():
    def __init__(self, fixed_params, fixed_params_vals, floated_params, floated_params_priors):
        self.fixed_params = fixed_params
        self.fixed_params_vals = fixed_params_vals
        self.floated_params = floated_params
        self.floated_params_priors = floated_params_priors

        self.all_params = ['N_bulge','N_disk','alpha','n','sigma','z0','beta_disk','beta_bulge']

        # Load the binned efficiency and data files
        self.PSR_data = np.load('../data/PSR_data.npy') + np.load('../data/PSR_data_3fgl.npy')
        self.omega_ijk = np.load('../data/omega_ijk.npy')


        self.setup_fixed_params()
        self.setup_prior_thetas()

    def setup_fixed_params(self):
        self.theta = np.zeros(len(self.all_params))
        theta_fixed = []
        for param in self.all_params:
            if param in self.fixed_params:
                param_pos = np.argwhere(np.array(self.fixed_params) == param)[0][0]
                param_all_pos = np.argwhere(np.array(self.all_params) == param)[0][0]
                theta_fixed += [self.fixed_params_vals[param_pos]]
                self.theta[param_all_pos] = self.fixed_params_vals[param_pos]

    def setup_prior_thetas(self):
        self.param_all_pos_ary = []
        self.theta_min = []
        self.theta_max = []
        for param in self.all_params:
            if param in self.floated_params:
                param_pos = np.argwhere(np.array(self.floated_params) == param)[0][0]
                param_all_pos = np.argwhere(np.array(self.all_params) == param)[0][0]
                self.param_all_pos_ary.append(param_all_pos)
                self.theta_min += [self.floated_params_priors[param_pos][0]]
                self.theta_max += [self.floated_params_priors[param_pos][1]]

        self.theta_interval = list(np.array(self.theta_max) - np.array(self.theta_min))

    def prior_cube(self, cube, ndim=1, nparams=1):
        """ Cube of priors - in the format required by MultiNest
        """
        for i in range(ndim):
            cube[i] = cube[i] * self.theta_interval[i] + self.theta_min[i]
        return cube

    def ll(self, theta, ndim = 1, nparams = 1):
        """ Log Likelihood in the format required by Multinest
        """
        theta_ll = self.theta
        for float_idx, float_item in enumerate(self.param_all_pos_ary):
            theta_ll[float_item] = theta[float_idx]
        return LL.ll(self.PSR_data, self.omega_ijk, *theta_ll)

    def perform_scan(self):
        # Run multinest
        n_params = len(self.floated_params)
        nlive = 200
        chains_dir = '/group/hepheno/smsharma/GCE-2FIG-Ben/run/chains/base'
        pymultinest_options = {'importance_nested_sampling': False,
                                'resume': False, 'verbose': True,
                                'sampling_efficiency': 'model',
                                'init_MPI': False, 'evidence_tolerance': 0.5,
                                'const_efficiency_mode': False}

        pymultinest.run(self.ll, self.prior_cube, n_params, outputfiles_basename = chains_dir, 
                        n_live_points = nlive, **pymultinest_options)

