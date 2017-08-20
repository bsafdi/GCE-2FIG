###############################################################################
# run.py
###############################################################################
#
# Perform scan to obtain best-fit disk and bulge PS population parameters
#
# Written: Siddharth Mishra-Sharma, Nick Rodd, and Ben Safdi 15 August 2017
# 
###############################################################################

import sys, os
sys.path.append("../python/")

import numpy as np
import pymultinest
from iminuit import Minuit
import corner

import likelihood
from minuit_functions import call_ll

class run_scan():
    def __init__(self, fixed_params, fixed_param_vals, floated_params, floated_param_priors, data_dir = '../data/',
                 Lmin = 1.0e31, Ns = 200, Nang = 1, smax_disk = 40, theta_mask=2.0, share_betas=False,
                 use_prior=False, Nang_prior=40, efficiency_long=False, efficiency_custom=None):
        """ Initialize scan class.

            :param fixed_params: Array of parameters to be held fixed
            :param fixed_param_vals: Array of values for parameters to be held fixed
            :param floated_params: Array of parameters to be floated
            :param floated_param_priors: Priors array for parameters to be floated
            :param data_dir: Directory containing the required maps
            :param Lmin: Minimum luminosity (erg s$^{-1}$)
            :param Lmax_disk: Maximum luminosity for disk (erg s$^{-1}$)
            :param Lmax_bulge: Maximum luminosity for bulge (erg s$^{-1}$)
            :param Ns: Number of integration point in z (kpc)
            :param Nang: Number of angular integration points
            :param smax_disk: How far to integrate out to for disk  (kpc)
            :param theta_mask: How many inner degrees to mask
            :param share_betas: Whether to float a single beta for disk and bulge
            :param use_prior: Whether to use prior on total number of sources
            :param Nang_prior: How many angular bins to use over full sky when using prior
            :param efficiency_long: Whether to use longitude-dependent efficiency
            :param efficiency_custom: A custom efficiency function
        """
        self.fixed_params = np.array(fixed_params)
        self.fixed_param_vals = np.array(fixed_param_vals)
        self.floated_params = np.array(floated_params)
        self.floated_param_priors = np.array(floated_param_priors)
        self.data_dir = data_dir

        self.Lmin = Lmin
        self.Ns = Ns
        self.Nang = Nang
        self.smax_disk = smax_disk
        self.theta_mask = theta_mask
        self.share_betas = share_betas

        self.efficiency_long = efficiency_long
        self.efficiency_custom = efficiency_custom

        self.use_prior = use_prior
        self.Nang_prior = Nang_prior
        
        self.all_params = np.array(['N_bulge','N_disk','alpha','n','sigma','z0','Lmax_disk','Lmax_bulge','beta_disk','beta_bulge'])
        
        if self.share_betas: # If using same beta for bulge and disk, remove the bulge param and use disk one in common
            
            self.all_params = self.all_params[self.all_params != 'beta_bulge']
            self.floated_params = self.floated_params[self.floated_params != 'beta_bulge']
            self.fixed_params = self.fixed_params[self.fixed_params != 'beta_bulge']

            self.all_params[self.all_params == 'beta_disk'] = 'beta'
            self.floated_params[self.floated_params == 'beta_disk'] = 'beta'

        self.load_data()
        self.setup_fixed_params()
        self.setup_prior_thetas()

    def load_data(self):
        """ Load the binned efficiency and data files
        """       
        self.PSR_data = np.load(self.data_dir + '/PSR_data.npy') + np.load(self.data_dir + '/PSR_data_3fgl.npy')
        
        if self.efficiency_custom is None: # If not using a custom efficiency
            if not self.efficiency_long: # If not using longitude-dependent efficiency
                omega_jk = np.load('../data/omega_jk.npy')
                self.omega_ijk = np.zeros((12, 12, 8))
                for i in range(12):
                    self.omega_ijk[i,:,:] = omega_jk
            else: # Load in longitude-dependent efficiency
                self.omega_ijk = np.load(self.data_dir + '/omega_ijk_2PC.npy')
        else: # Load in custom efficiency if this is specified
            self.omega_ijk = self.efficiency_custom

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
        """ Cube of priors in the format required by MultiNest
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

        if self.share_betas and ('beta' in self.floated_params):
            theta_ll = np.append(theta_ll, theta_ll[-1])  

        ll_val =  likelihood.ll(self.PSR_data, self.omega_ijk,                  
                        *theta_ll, 
                        Lmin = self.Lmin, 
                        Ns = self.Ns, Nang = self.Nang, smax_disk = self.smax_disk, theta_mask = self.theta_mask,
                        use_prior = self.use_prior, Nang_prior = self.Nang_prior)

        return ll_val

    def perform_scan_multinest(self, nlive=100, chains_dir='/group/hepheno/smsharma/GCE-2FIG-bsafdi/run/chains/analysis/', resume=False):
        """ Perform a scan with MultiNest
        """
        self.make_dirs([chains_dir])
        n_params = len(self.floated_params)
        pymultinest_options = {'importance_nested_sampling': False,
                                'resume': resume, 'verbose': True,
                                'sampling_efficiency': 'model',
                                'init_MPI': False, 'evidence_tolerance': 0.5,
                                'const_efficiency_mode': False}

        pymultinest.run(self.ll, self.prior_cube, n_params, outputfiles_basename = chains_dir, 
                        n_live_points = nlive, **pymultinest_options)

    def perform_scan_minuit(self, verbose=1):
        """ Perform a scan with minuit
        """
        limit_dict = {}
        init_val_dict = {}
        step_size_dict = {}
        for iparam, param in enumerate(self.floated_params):
            limit_dict['limit_'+param] = (self.floated_param_priors[iparam][0],self.floated_param_priors[iparam][1])
            init_val_dict[param] = (self.floated_param_priors[iparam][0] + self.floated_param_priors[iparam][1])/2.
            step_size_dict['error_'+param] = 1.0
        other_kwargs = {'print_level': verbose, 'errordef': 1}
        z = limit_dict.copy()
        z.update(other_kwargs)
        z.update(limit_dict)
        z.update(init_val_dict)
        z.update(step_size_dict)
        f = call_ll(len(self.floated_params),self.ll,self.floated_params)
        self.m = Minuit(f,**z)
        self.m.migrad()
        self.max_LL = - self.m.fval

    def plot_corner(self, chains_dir, labels):
        # Load the samples
        chain_file = chains_dir + 'post_equal_weights.dat'
        self.samples = np.array(np.loadtxt(chain_file)[:, :-1])

        # Now make a triangle plot using corner
        corner.corner(self.samples, labels=labels, smooth=1.5,
                      smooth1d=1, quantiles=[0.16, 0.5, 0.84], show_titles=True,
                      title_fmt='.2f', title_args={'fontsize': 14},
                      range=[1 for _ in range(len(self.floated_params))],
                      plot_datapoints=False, verbose=False)

    def get_lge(self, chains_dir):
        # Load the samples
        a = pymultinest.Analyzer(n_params=len(self.floated_params),
                                     outputfiles_basename=chains_dir)
        s = a.get_stats()

        lge = s['nested sampling global log-evidence']
        lge_err = s['nested sampling global log-evidence error']

        return lge, lge_err

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
