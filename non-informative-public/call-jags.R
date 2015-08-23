



SEED <- 0

crop <- FALSE # If desired
star <- 0 #!!


# MCMC settings
#ni <- 100; nb <- 20; nt <- 5; nc <- 1
#ni <- 1000; nb <- 20; nt <- 5; nc <- 1
ni <- 50000; nb <- 25000; nt <- 20; nc <- 1



#load necessary packages
library("lme4")
library("bayesm")
library("R2jags")
library("splines")
library("dplyr")





#######################################
#######################################
#Load & filter data

pt.data<-read.csv("simulation-data/pt-data-sim.csv")
psa.data.all<-read.csv("simulation-data/psa-data-sim.csv")
data.use.all<-read.csv("simulation-data/bx-data-sim.csv")
#this contains one record per annual interval for each patient until surgery or censoring
#data frames with '.all' suffix will be cropped in some cases to for "leave one out" JAGS fits.

#Crop out data for "leave-one-out" model fits.
data.use.star<- filter(data.use.all, subj==star)
bx.inds <- which(data.use.star$bx.here==1)
last.bx.ind <- bx.inds[length(bx.inds)-1] #Take all data since the 2nd to last biopsy as "new data."
last.age <- data.use.star$age[last.bx.ind]


if(!crop | length(last.age)==0 ) last.age <- 0 #there could be some PSA measurements before you're enrolled in this biopsy study.
psa.data <- filter(psa.data.all, !(subj==star & age>=last.age))
data.use <- filter(data.use.all, !(subj==star & age>=last.age))




(n<-dim(pt.data)[1]) #there are 1000 patients
#This matrix will always be fully intact, but psa and bx data may not be.




#######################################
#######################################
#Generate covariates and datavectors to be sent to JAGS



#function to get natrual spline basis
#to be used if IOP_BX or IOP_RRP is true
get.ns.basis<-function(obs.data,knots){
#	knots<-quantile(obs.data,p=c(0.25,0.5,0.75))
	od.k1<- obs.data-knots[1]
	od.k1[od.k1<0]<-0
	od.k2<- obs.data-knots[2]
	od.k2[od.k2<0]<-0
	od.k3<- obs.data-knots[3]
	od.k3[od.k3<0]<-0
	return(as.vector((od.k1^3 - od.k3^3)/(knots[3]-knots[1]) - (od.k2^3 - od.k3^3)/(knots[3]-knots[2])))}




#get observed latent class
eta.data<-pt.data$obs.eta
table(eta.data) #107 in each
(n_eta_known<-sum(!is.na(eta.data))) #214


#########
# PSA model

(n_obs_psa<-dim(psa.data)[1])
Y<-psa.data$log.psa
summary(Y)
subj_psa<-psa.data$subj

#covariates with random effects
Z.data<-as.matrix(cbind(rep(1,n_obs_psa), psa.data$age.std))
(d.Z<-dim(Z.data)[2])
round(apply(Z.data,2,summary),2)

#covariates with only fixed effects
X.data<-as.matrix(cbind(psa.data$std.vol))
(d.X<-dim(X.data)[2])
summary(X.data)
#########


#########
# Reclassification (RC) logistic regression outcome model

rc.data<-data.use[data.use$bx.here==1 & !is.na(data.use$bx.here),] #only use records where a biopsy occurred
(n_rc<-dim(rc.data)[1])
RC<-as.numeric(rc.data$rc)
subj_rc<-rc.data$subj

#covariates influencing risk of reclassification
W.RC.data<-as.matrix(cbind(rep(1,n_rc),  rc.data$age.std, rc.data$time, rc.data$time.ns, rc.data$sec.time.std ))
(d.W.RC<-dim(W.RC.data)[2])
round(apply(W.RC.data,2,summary) ,2)
#########



#########
#biopsy (BX) logistic regression observation model

#If these are NULL, they effectively won't be in the list passed to JAGS
bx.data <-
n_bx <-
BX <-
subj_bx<-
W.BX.data<-
d.W.BX <- NULL

if(IOP_BX){
	bx.data<-data.use[!is.na(data.use$bx.here),] # remove patients who have already had RC observed but haven't had surgery or been censored
	(n_bx<-dim(bx.data)[1])
	BX<-as.numeric(bx.data$bx.here) #indicator of bx
	subj_bx<-bx.data$subj


	W.BX.data<-as.matrix(cbind(rep(1,n_bx), bx.data$age.std, bx.data$age.ns, ns(bx.data$time,4), bx.data$num.prev.bx, ns(bx.data$sec.time.std,4)  ))
	(d.W.BX<-dim(W.BX.data)[2]) 
	round(apply(W.BX.data,2,summary),2)
}
#########




#########
# Surgery (RRP) logistic regression observation model
#logistic regression for RRP
#this uses all records, because patients always at risk of choosing surgery

RRP<-
n_rrp<-
subj_rrp<-
W.RRP.data<-
d.W.RRP<- NULL

if(IOP_RRP){
	RRP<-as.numeric(data.use$rrp)
	(n_rrp<-dim(data.use)[1])
	subj_rrp<-data.use$subj


	W.RRP.data<-as.matrix(cbind(rep(1,n_rrp), data.use$age.std, data.use$age.ns, ns(data.use$time,4), ns(data.use$sec.time.std,3) , data.use$num.prev.bx.rrp, data.use$prev.G7)) #
	(d.W.RRP<-dim(W.RRP.data)[2])
	round(apply(W.RRP.data,2,summary) ,2)
}
#########





#######################################
#######################################
# Initialize parameters for JAGS
# 	(and other JAGS helper functions)

#lmer fit for initializing variance parameter in JAGS
mod.lmer<-lmer(log.psa~ std.vol + (1+ age.std |id), data=psa.data)
(var_vec <- apply(coef(mod.lmer)$id, 2, var)[1:d.Z])
(var_vec <- c(var_vec[2], var_vec[1]))

#bundle data for call to JAGS
#this is observed data and constant variables that have already been assigned values (e.g. number of class K=2, number of subjects n, etc.)
K<-2

jags_data<-list(K=K, n=n, eta.data=eta.data, n_eta_known=n_eta_known,
	Y=Y, subj_psa=subj_psa, Z=Z.data, X=X.data, d.Z=d.Z, d.X=d.X, I_d.Z=diag(d.Z), n_obs_psa=n_obs_psa,
	RC=RC, n_rc=n_rc, subj_rc=subj_rc, W.RC=W.RC.data, d.W.RC=d.W.RC,
	IOP_BX=IOP_BX, IOP_RRP=IOP_RRP)

if(IOP_BX) jags_data <- c(jags_data,
		list(BX=BX, n_bx=n_bx, subj_bx=subj_bx, W.BX=W.BX.data, d.W.BX=d.W.BX))
if(IOP_RRP) jags_data <- c(jags_data, 
		list(RRP=RRP, n_rrp=n_rrp, subj_rrp=subj_rrp, W.RRP=W.RRP.data, d.W.RRP=d.W.RRP))


#initialize model
#this is to set initial values of parameters
#note that not all "parameters" need to be initialized. specifically, don't initialize random effects, but do need to set initial values for mean and covariance of random effects
# also need to set initial values for latent variables that are not observed (here, eta)

inits <- function(){
		
	p_eta<-rbeta(1,1,1)

	eta.hat<-pt.data$rc[is.na(eta.data)]

	xi<-c(min(rlnorm(1),100), min(rlnorm(1),100))
	mu_raw<-as.matrix(cbind(rnorm(d.Z),rnorm(d.Z)))
	Tau_B_raw<-rwishart((d.Z+1),diag(d.Z)*var_vec)$W
	sigma_res<-min(rlnorm(1),1)

	beta<-rnorm(d.X)

	gamma.RC<-rnorm((d.W.RC+1), mean=0, sd=0.1) #last coefficient is effect of eta=1

	gamma.BX <- 
	gamma.RRP <- NULL
	if(IOP_BX) gamma.BX<-rnorm((d.W.BX+1), mean=0, sd=0.1) #last coefficient is effect of eta=1
	if(IOP_RRP) gamma.RRP<-c(rnorm((d.W.RRP+2), mean=0, sd=0.01))  #here, include interaction with last prediction and eta=1

	out<-list(p_eta=p_eta, eta.hat=eta.hat, xi=xi, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, sigma_res=sigma_res, beta=beta,
		gamma.RC=gamma.RC)
	if(IOP_BX) out$gamma.BX <- gamma.BX
	if(IOP_RRP) out$gamma.RRP <- gamma.RRP
	
	out

}



# parameters to track
params <- c("p_eta", "eta.hat", "mu_int", "mu_slope", "sigma_int", "sigma_slope", "sigma_res", "rho_int_slope", "cov_int_slope", "b.vec", "beta",
	"gamma.RC",
	"p_rc")
if(IOP_BX) params <- c(params, "gamma.BX", "p_bx")
if(IOP_RRP) params <- c(params, "gamma.RRP", "p_rrp")	
#you may not need to monitor p_bx, p_rc, and p_rrp. taking them out of the list should improve computing time a bit

# MCMC settings
#ni, nb, nt, and nc are now set in separate files.


model.file.IOP <- paste0(
	'model-for-jags-',
	c('N')[!IOP_BX],'IOP_BX-',
	c('N')[!IOP_RRP],'IOP_RRP',
	'.txt'
	)

do.one<-function(seed){
	set.seed(seed)	
	outj<-jags(jags_data,
		inits=inits,
		parameters.to.save=params,
		model.file=model.file.IOP,
		n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

	return(outj$BUGSoutput)
}



outJAGS <- do.one(seed=SEED) #final output of this script


#Reduce filesize of results (make new object, outJAGS2)
#Parameters to ignore
params2ignore<-c('b.vec','rho_int_slope','eta.hat','p_rc')
if(IOP_BX) params2ignore <- c(params2ignore, 'p_bx')
if(IOP_RRP) params2ignore <- c(params2ignore, 'p_rrp')


outJAGS2_sl <- outJAGS$sims.list[!names(outJAGS$sims.list) %in% params2ignore]
outJAGS2 <- outJAGS[!names(outJAGS) %in% c('sims.matrix','sims.array')]
outJAGS2_sl$eta_hat_means <- apply(outJAGS$sims.list$eta.hat,2,mean)
if(crop){ #if `crop`, save subj-specific posterior for subj star. otherwise don't.
	outJAGS2_sl$b.vec.star <- outJAGS$sims.list$b.vec[,star,]
	outJAGS2_sl$eta.hat.star <- NULL
	if(star > 214)
		outJAGS2_sl$eta.hat.star <- outJAGS$sims.list$eta.hat[,star - 214]
}

outJAGS2$sims.list <- outJAGS2_sl


len.sim<-length(outJAGS2$sims.list$p_eta)
saveRDS(outJAGS2,file=paste0(save_path,Sys.Date(),'_posterior_IOP_BX-',IOP_BX,'_IOP_RRP-',IOP_RRP,'_SEED-',SEED,'_star-',star,'_crop-',crop,'_nsim-',len.sim,'.rds'))
str(outJAGS2$sims.list)
summary(outJAGS2$sims.list$mu_slope)









