# Code for submitting jobs and running analysis on a Sun Grid Engine cluster

# To begin, type 'sh run_pipeline.sh' into the command line.

# This code can also be run in parallel to increase the
# size of the posterior draw. However, for simplicity of
# demonstration and easy of reproducibility across computing
# systems, we only present a serial example (without parallelization).

#Step 1: Run MCMC
qsub -N runMCMC -V -l mf=4G,h_vmem=4G -cwd -b y 'R CMD BATCH --no-save call_jags.R'

#Step 2: Run Importance Sampling
qsub -N runIS -hold_jid runMCMC -V -l mf=4G,h_vmem=4G -cwd -b y 'R CMD BATCH --no-save fit_IS.R fit_IS.Rout'


