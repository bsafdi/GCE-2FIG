import os, sys
import numpy as np

Nb_ary = 10**np.linspace(0,7,20)

for Nb in Nb_ary:

    batch1='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=12
#SBATCH -t 48:00:00
#SBATCH --mem=36GB
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=smsharma@princeton.edu
#SBATCH -p dept

source activate venv_py27

export MKL_LP64_ILP64=ilp64
source /opt/intel/compilers_and_libraries_2016.2.181/linux/bin/compilervars.sh intel64
source /opt/intel/impi/5.0.3.048/bin64/mpivars.sh

export LDSHARED="icc -shared" CC=icc

# module load openmpi/gcc/1.6.5/64

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/group/hepheno/heptools/MultiNest/lib/

cd /group/hepheno/smsharma/GCE-2FIG-bsafdi/run/
'''
    batch2 ='Nb='+str(Nb)+'\n'
    batch3 = '''
mpiexec.hydra -n 12 python run_nd_nb_z0_beta_prior_profile.py --Nb $Nb\n
'''

    batchn = batch1+batch2+batch3
    fname = "./batch/run_Nb_profile-"+str(Nb)+".batch"
    f=open(fname, "w")
    f.write(batchn)
    f.close()
    os.system("sbatch "+fname);
