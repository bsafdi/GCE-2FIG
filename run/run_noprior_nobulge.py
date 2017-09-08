from run import run_scan

fixed_params = ['n','sigma','alpha','Lmax_disk', 'Lmax_bulge','N_bulge']
fixed_param_vals = [2.35,1.528,2.6,1.e36,1.e36,0.]

floated_params = ['N_disk', 'z0', 'beta']
floated_param_priors = [[0.,1.e7],[0.01,2.],[1.1,3.]]

rs_nd = run_scan(fixed_params, fixed_param_vals, floated_params, floated_param_priors, Ns = 200, Nang = 5, share_betas=True, use_prior=False)
rs_nd.perform_scan_multinest(nlive=500, chains_dir='chains/np_nb_nlive500/')
