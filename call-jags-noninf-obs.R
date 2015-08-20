rm(list=ls())
# setwd("/Users/ryc/Dropbox/inhealth/prediction-model")
# setwd("/Users/ryc/GitHub/prostate_surveillance")



# !! When you do the importance sampling comparison for *new data* for person with existing data, don't just do last PSA, go back to last biopsy, and treat all data from that day on (including the biopsy) as "newly acquired" data.
# !! Use age field (not age.std) to separate newly acquired from previously acquired data.
# Make a plot of "new patient" and plot of "new data on existing patient".
# Run on actual data eventually.


#load necessary packages
library("lme4")
library("bayesm")
library("R2jags")
library("splines")
library("dplyr")


#get data
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



#function to get natrual spline basis
get.ns.basis<-function(obs.data,knots){
#	knots<-quantile(obs.data,p=c(0.25,0.5,0.75))
	od.k1<- obs.data-knots[1]
	od.k1[od.k1<0]<-0
	od.k2<- obs.data-knots[2]
	od.k2[od.k2<0]<-0
	od.k3<- obs.data-knots[3]
	od.k3[od.k3<0]<-0
	return(as.vector((od.k1^3 - od.k3^3)/(knots[3]-knots[1]) - (od.k2^3 - od.k3^3)/(knots[3]-knots[2])))}


#Before call to JAGS, get the data into simple matrices and vectors to send to JAGS

(n<-dim(pt.data)[1]) #there are 1000 patients
#This matrix will always be fully intact, but psa and bx data may not be.

#get observed latent class
eta.data<-pt.data$obs.eta
table(eta.data) #107 in each
(n_eta_known<-sum(!is.na(eta.data))) #214


#Before call to JAGS, get the data into simple matrices and vectors to send to JAGS

#latent class model
#no regression model here

#PSA model
(n_obs_psa<-dim(psa.data)[1])
	#this is the number of PSA observations we have, so it is >>number subjects
	#I will loop through all 1:n_obs_psa observations in the JAGS code

Y<-psa.data$log.psa
summary(Y)
subj_psa<-psa.data$subj #unique id that corresponds to the vector eta.data and will be used to index the random effects


#covariates with random effects
Z.data<-as.matrix(cbind(rep(1,n_obs_psa), psa.data$age.std))

#this is the design matrix for the random effects
#age.std is standardized age, i.e., centered at mean of all ages for PSA observations and divided by the sd of those ages
#psa.age.basis is the 3rd basis fcn, h_3(), in the write-up. it corresponds to the "curvature"

(d.Z<-dim(Z.data)[2])
round(apply(Z.data,2,summary),2) #this is just here to check that I defined things correctly, have pulled in the right data


#covariates with only fixed effects
X.data<-as.matrix(cbind(psa.data$std.vol))
(d.X<-dim(X.data)[2])
summary(X.data)


#outcome model (logistic regression for reclassification)
rc.data<-data.use[data.use$bx.here==1 & !is.na(data.use$bx.here),]

(n_rc<-dim(rc.data)[1])
RC<-as.numeric(rc.data$rc)
#table(RC) #300

subj_rc<-data.use$subj

#covariates influencing risk of reclassification
W.RC.data<-as.matrix(cbind(rep(1,n_rc),  rc.data$age.std, rc.data$time, rc.data$time.ns, rc.data$sec.time.std)) #this last predictor is a measure of secular time (biopsy grading trends changed over time)   
(d.W.RC<-dim(W.RC.data)[2])
round(apply(W.RC.data,2,summary),2)


##get starting values, other functions necessary for call to JAGS


#lmer fit for initializing parameters
#do this to get the starting value for a variance paramter in JAGS
mod.lmer<-lmer(log.psa~ std.vol + (1+ age.std |id), data=psa.data)
(var_vec <- apply(coef(mod.lmer)$id, 2, var)[1:d.Z]) #not sure why these aren't printed in the right order...
(var_vec<- c(var_vec[2], var_vec[1])) #fixing order

#bundle data for call to JAGS
#this is observed data and constant variables that have already been assigned values (e.g. number of class K=2, number of subjects n, etc.)
K<-2

relabel_consecutive<-function(x){
	fx<-as.factor(x)
	rlfx<-mapvalues(fx,
		from=levels(fx),
		to=1:length(levels(fx)))
	as.numeric(rlfx)
}
subj_psa_consecutive <- relabel_consecutive(subj_psa)
subj_rc_consecutive <- relabel_consecutive(subj_rc)

jags_data<-list(K=K,
	n=n,
	eta.data=eta.data,
	n_eta_known=n_eta_known,
	n_obs_psa=n_obs_psa,
	Y=Y,
	subj_psa=subj_psa_consecutive, #!! new
	Z=Z.data,
	X=X.data,
	d.Z=d.Z,
	d.X=d.X,
	I_d.Z=diag(d.Z),
	n_rc=n_rc,
	RC=RC,
	subj_rc=subj_rc_consecutive,#!! new
	W.RC=W.RC.data,
	d.W.RC=d.W.RC)


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

	gamma.BX<-rnorm((d.W.BX+1), mean=0, sd=0.1) #last coefficient is effect of eta=1
	gamma.RC<-rnorm((d.W.RC+1), mean=0, sd=0.1) #ditto
	gamma.RRP<-c(rnorm((d.W.RRP+2), mean=0, sd=0.01))  #here, include interaction with last prediction and eta=1

	list(p_eta=p_eta, eta.hat=eta.hat, xi=xi, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, sigma_res=sigma_res, beta=beta, gamma.BX=gamma.BX, gamma.RC=gamma.RC, gamma.RRP=gamma.RRP)
}



# parameters to track
params <- c("p_eta", "eta.hat", "mu_int", "mu_slope", "sigma_int", "sigma_slope", "sigma_res", "rho_int_slope", "cov_int_slope", "b.vec", "beta", "gamma.RC")

# MCMC settings
#ni, nb, nt, and nc are now set in separate files.


do.one<-function(seed){
	set.seed(seed)	
	outj<-jags(jags_data, inits=inits, parameters.to.save=params, model.file="model-for-jags-inf-obs.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

	return(outj$BUGSoutput)
}










