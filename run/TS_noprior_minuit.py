###############################################################################
# TS_noprior_minuit.py
###############################################################################
# 
# Evaluate the TS for the bulge using Minuit without the prior
#
###############################################################################

# Import run module
from run import run_scan


# Run the scan with the bulge
fixed_params = ['n','sigma','alpha','Lmax_disk','Lmax_bulge']
fixed_param_vals = [2.35,1.528,2.6,1.e36,1.e36]
floated_params = ['N_bulge','N_disk','z0','beta']
floated_param_priors = [[0.,1.e6],[0.,1.e6],[0.01,2.],[1.1,3.]]

rs_wb = run_scan(fixed_params, fixed_param_vals, floated_params, 
                 floated_param_priors, Nang=1, share_betas=True, 
                 use_prior=False)
rs_wb.perform_scan_minuit(verbose=0)
wb_ll = rs_wb.max_LL


# Run the scane without the bulge
fixed_params = ['n','sigma','alpha','Lmax_disk','Lmax_bulge','N_bulge']
fixed_param_vals = [2.35,1.528,2.6,1.e36,1.e36,0.]
floated_params = ['N_disk','z0','beta']
floated_param_priors = [[0.,1.e6],[0.01,2.],[1.1,3.]]

rs_nb = run_scan(fixed_params, fixed_param_vals, floated_params,
                 floated_param_priors, Nang=1, share_betas=True,
                 use_prior=False)
rs_nb.perform_scan_minuit(verbose=0)
nb_ll = rs_nb.max_LL


# Print the TS
print 2.*(wb_ll-nb_ll)
