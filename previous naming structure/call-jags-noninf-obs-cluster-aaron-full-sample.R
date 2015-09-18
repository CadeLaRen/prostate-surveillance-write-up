# cd /home/bst/student/afisher/inHealth_prostate
# qsub -N fullJAGS -t 1-200 -V -l mf=5G,h_vmem=5G -cwd -b y R CMD BATCH --no-save call-jags-noninf-obs-cluster-aaron-full-sample.R
# qsub -N leaveOutJAGS -t 215-1000 -V -l mf=5G,h_vmem=5G -cwd -b y R CMD BATCH --no-save call-jags-noninf-obs.R
# rm leaveOutJAGS.e*
# rm leaveOutJAGS.o*





#import environment variable, used for running multiple chains in parallel
(SEED<-as.numeric(Sys.getenv("SGE_TASK_ID")))
if(is.na(SEED)) SEED <- 0

star <- 0

# MCMC settings
#ni <- 100; nb <- 20; nt <- 5; nc <- 1
#ni <- 1000; nb <- 20; nt <- 5; nc <- 1
ni <- 50000; nb <- 25000; nt <- 20; nc <- 1

cat('',file=paste0('batches/noninf-obs/1/',Sys.Date(),'_JAGS_started_SEED_',SEED,'_star_',star,'.txt'))

source('call-jags-noninf-obs.R') #returns object `outJAGS`

len.sim<-length(outJAGS$sims.list$p_eta)
saveRDS(outJAGS,file=paste0('batches/noninf-obs/1/',Sys.Date(),'_posterior_noninf-obs_SEED-',SEED,'_star-',star,'_nsim-',len.sim,'.rds'))
str(outJAGS$sims.list)
summary(outJAGS$sims.list$mu_slope)

