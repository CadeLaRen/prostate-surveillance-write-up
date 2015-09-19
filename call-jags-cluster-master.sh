# cd /home/bst/student/afisher/inHealth_prostate

#args to R correspond to: `IOP_BX`, `IOP_SURG`, and `leave_one_out` respectively


#Fit on entire sample
qsub -N IOPleave1out -t 1-200 -V -l mf=4G,h_vmem=4G -cwd -b y 'R CMD BATCH --no-save "--args TRUE TRUE TRUE" call-jags-cluster-master.R IOPleave1out.Rout'
qsub -N IOPall -t 1-200 -V -l mf=4G,h_vmem=4G -cwd -b y 'R CMD BATCH --no-save "--args TRUE TRUE FALSE" call-jags-cluster-master.R IOPall.Rout'


qstat | wc -l
qstat

rm IOPall.o*
rm IOPall.e*

rm leave1out.o*
rm leave1out.e*

tail IOPall.Rout -n 100
tail IOPleave1out.Rout -n 100

# For checking output from other interactive sessions
	# tail call-jags-cluster-master.R -n 30








###################################
###################################
###################################
# Versions for non-informative observation models

# qsub -N noBX_yesRRP -t 1-50 -V -l mf=4G,h_vmem=4G -cwd -b y 'R CMD BATCH --no-save "--args FALSE TRUE FALSE" call-jags-cluster-master.R noBX_yesRRP.Rout'
# qsub -N yesBX_noRRP -t 1-50 -V -l mf=4G,h_vmem=4G -cwd -b y 'R CMD BATCH --no-save "--args TRUE FALSE FALSE" call-jags-cluster-master.R yesBX_noRRP.Rout'
# qsub -N noBX_noRRP -t 1-50 -V -l mf=4G,h_vmem=4G -cwd -b y 'R CMD BATCH --no-save "--args FALSE FALSE FALSE" call-jags-cluster-master.R noBX_noRRP.Rout'

# rm noBX_yesRRP.e*
# rm noBX_yesRRP.o*
# rm yesBX_noRRP.o*
# rm yesBX_noRRP.e*
# rm noBX_noRRP.e*
# rm noBX_noRRP.o*

# tail noBX_yesRRP.Rout -n 100
# tail yesBX_noRRP.Rout -n 100
# tail noBX_noRRP.Rout -n 100
