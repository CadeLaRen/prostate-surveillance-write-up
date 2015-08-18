# cd /home/bst/student/afisher/inHealth_prostate
# #qsub -N fullJAGS -t 1-200 -V -l mf=5G,h_vmem=5G -cwd -b y R CMD BATCH --no-save call-jags-inf-obs-cluster-aaron-full-sample.R
# qsub -N oneVisitOut -t 215-1000 -V -l mf=5G,h_vmem=5G -cwd -b y R CMD BATCH --no-save call-jags-inf-obs-cluster-aaron-leave-1visit-out.R
# rm oneVisitOut.o*
# rm oneVisitOut.e*

# tail call-jags-inf-obs-cluster-aaron-leave-1visit-out.R -n 30
# tail call-jags-inf-obs-cluster-aaron-leave-1visit-out.Rout



#import environment variable, used for running multiple chains in parallel
(SEED<-as.numeric(Sys.getenv("SGE_TASK_ID")))
if(is.na(SEED)) SEED <- 0

crop <- TRUE
star <- SEED

# MCMC settings
#ni <- 100; nb <- 20; nt <- 5; nc <- 1
#ni <- 1000; nb <- 20; nt <- 5; nc <- 1
ni <- 50000; nb <- 25000; nt <- 20; nc <- 1

cat('',file=paste0('leaveOneOut/',Sys.Date(),'_started_inf-obs_SEED-',SEED,'_star-',star,'_crop-',crop,'.txt'))

source('call-jags-inf-obs.R') #returns function `do.one`
outJAGS <- do.one(seed=SEED) #final output of this script


#Reduce filesize of results (make new object, outJAGS2)
#Parameters to ignore
params2ignore<-c('b.vec','rho_int_slope','eta.hat','p_bx','p_rc','p_rrp')


outJAGS2_sl <- outJAGS$sims.list[!names(outJAGS$sims.list) %in% params2ignore]
outJAGS2 <- outJAGS[!names(outJAGS) %in% c('sims.matrix','sims.array')]
outJAGS2_sl$eta_hat_means <- apply(outJAGS$sims.list$eta.hat,2,mean)
if(crop){
	outJAGS2_sl$b.vec.star <- outJAGS$sims.list$b.vec[,star,]
	outJAGS2_sl$eta.hat.star <- NULL
	if(star > 214)
		outJAGS2_sl$eta.hat.star <- outJAGS$sims.list$eta.hat[,star - 214]
}

outJAGS2$sims.list <- outJAGS2_sl


len.sim<-length(outJAGS2$sims.list$p_eta)
saveRDS(outJAGS2,file=paste0('leaveOneOut/',Sys.Date(),'_posterior_inf-obs_SEED-',SEED,'_star-',star,'_crop-',crop,'_nsim-',len.sim,'.rds'))
str(outJAGS2$sims.list)
summary(outJAGS2$sims.list$mu_slope)

