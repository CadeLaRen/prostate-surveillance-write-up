#setwd("/Users/ryc/Dropbox/inhealth/prediction-model-final/sim-data")
#setwd("/Users/ryc/GitHub/prostate_surveillance")

#import environment variable, used for running multiple chains in parallel
(SEED<-as.numeric(Sys.getenv("SGE_TASK_ID")))
if(is.na(SEED)) SEED <- 0

crop <- FALSE
star <- SEED

# MCMC settings
#ni <- 100; nb <- 20; nt <- 5; nc <- 1
#ni <- 1000; nb <- 20; nt <- 5; nc <- 1
ni <- 50000; nb <- 25000; nt <- 20; nc <- 1

source('call-jags-inf-obs.R') #returns function `do.one`
outJAGS <- do.one(seed=SEED) #final output of this script

for(j in 1:length(outJAGS$sims.list)){
	write.csv(outJAGS$sims.list[[j]], paste("jags-prediction-inf-obs-", names(out$sims.list)[j],"-",seed,".csv",sep=""))
}

