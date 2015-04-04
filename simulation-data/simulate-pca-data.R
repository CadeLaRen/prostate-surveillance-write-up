rm(list=ls())
setwd("/Users/ryc/Dropbox/inhealth/prediction-model/sim-data")


library(MASS)

set.seed(1)

#function to get natural spline basis
expit<-function(x){return(exp(x)/(1+exp(x)))}

get.ns.basis<-function(obs.data,knots){
	od.k1<- obs.data-knots[1]
	od.k1[od.k1<0]<-0
	od.k2<- obs.data-knots[2]
	od.k2[od.k2<0]<-0
	od.k3<- obs.data-knots[3]
	od.k3[od.k3<0]<-0
	return(as.vector((od.k1^3 - od.k3^3)/(knots[3]-knots[1]) - (od.k2^3 - od.k3^3)/(knots[3]-knots[2])))}


#parameter values to use for estimation (similar to posterior estimates)
p_eta <- 0.4

mu_int <- c(1.3, 1.6)
mu_slope <- c(0.3,0.5)
mu_spline <- c(-0.11, -0.07)
mu.mat <- as.matrix(rbind(mu_int, mu_slope, mu_spline))

Sigma <- matrix(c(0.65^2, 0.2, -0.2, 0.2, 0.5^2, -0.1, -0.2, -0.1, 0.4^2), nrow=3, ncol=3)
sigma_res <- 0.27

beta <- c(0.4,-0.3)

gamma <- c(-3.7, 0.5, 0.02, 1.8)



#
n <- 1000
id <- c(1:1000)

ages.dx <- rnorm(n, mean=65.5, sd=5.5) #from data

pt.data<-as.data.frame(cbind(id,ages.dx))
names(pt.data) <- c("id","age.dx")

#latent class
pt.data$eta.true <- eta.true <- rbinom(n,1,p_eta)

pt.data$obs.eta <- obs.eta <- rbinom(n,1,0.2) #no information in observation
(n_eta_known <- sum(obs.eta)) #229

write.csv(pt.data,"pt-data-sim.csv")

eta.data <- vector(length=n)
for(i in 1:n){
	if(obs.eta[i]==1){eta.data[i]<-eta.true[i]}
	else{eta.data[i]<-NA}	}

	
ordered <- order(eta.data)
eta.data <- eta.data[ordered]
pt.ordered <- pt.data[ordered,]

ids.ordered <- unique(pt.ordered$id)

write.csv(eta.data,"eta-data-sim.csv")
write.csv(pt.ordered, "pt-ordered-sim.csv")

#reclassification data
n_bx_ps <- 5 #number of biopsies per subject ("ps"), synonomous to number of years of follow-up 
bx.data <- as.data.frame(cbind(rep(1:n,n_bx_ps)))
names(bx.data) <- "id"

bx.time <- rep(1,n)
for(index in 2:n_bx_ps){bx.time <- c(bx.time, rep(index,n))}


bx.data$bx.time <- bx.time + runif((n*n_bx_ps),-0.5,0.5)
bx.data$bx.age <- vector(length=(n*n_bx_ps))
for(j in 1:(n*n_bx_ps)){
	bx.data$bx.age[j] <- bx.data$bx.time[j] + pt.data$age.dx[pt.data$id==bx.data$id[j]]}
bx.data$std.age <- (bx.data$bx.age-mean(bx.data$bx.age))/sd(bx.data$bx.age)


bx.data$rc <- vector(length=(n*n_bx_ps))
pt.data$rc <- pt.ordered$rc <- rep(0,n)
for(j in 1:(n*n_bx_ps)){
	cov.vec <- c(1, bx.data$std.age[j], bx.data$bx.time[j], pt.data$eta.true[pt.data$id==bx.data$id[j]])
	bx.data$rc[j] <- rbinom(1,1, expit(sum(cov.vec*gamma)))
	if(bx.data$rc[j]==1){pt.data$rc[pt.data$id==bx.data$id[j]]<-1
		pt.ordered$rc[pt.ordered$id==bx.data$id[j]]<-1}}


bx.data$rm <- rep(0,(n*n_bx_ps))
for(i in 1:n){
	time.rc<-time.cens<-NULL 
	if(pt.data$rc[i]==1){
		time.rc <- min(bx.data$bx.time[bx.data$rc==1 & bx.data$id==pt.data$id[i]])
		bx.data$rm[bx.data$id==pt.data$id[i] & bx.data$bx.time>time.rc] <- 1 } 
		else{
			time.cens <- c(1:5) %*% rmultinom(1,1,p=c(0.05,0.2,0.4,0.2,0.15))
			if(time.cens<5){
				bx.data$rm[bx.data$id==pt.data$id[i] & bx.data$bx.time > (time.cens[1,1]+0.5)] <- 1} }}


bx.data<-bx.data[bx.data$rm==0,]
bx.data$rm<-NULL

(n_obs_bx <- dim(bx.data)[1]) #3055
table(bx.data$rc) #290

bx.data$subj<-vector(length=n_obs_bx)
for(j in 1:n_obs_bx){
	bx.data$subj[j] <- c(1:n)[ids.ordered==bx.data$id[j]]}


write.csv(bx.data,"bx-data-sim.csv")
write.csv(pt.data,"pt-data-sim.csv")
write.csv(pt.ordered, "pt-ordered-sim.csv")



#psa data

psa.time <- seq(-1, max(bx.data$bx.time[bx.data$id==1]), 0.25)
psa.id <- rep(1, length(seq(-1, max(bx.data$bx.time[bx.data$id==1]), 0.25)))
for(i in 2:n){psa.time <- c(psa.time, seq(-1, max(bx.data$bx.time[bx.data$id==i]), 0.25) )
	psa.id <- c(psa.id, rep(i, length(seq(-1, max(bx.data$bx.time[bx.data$id==i]), 0.25)) ) ) }


psa.data <- as.data.frame(cbind(psa.id, psa.time))
names(psa.data) <- c("id", "psa.time")
(n_obs_psa<-dim(psa.data)[1]) #16803

psa.data$psa.time <- psa.data$psa.time + runif(n_obs_psa, min=-0.125, max=0.125)

psa.data$psa.age <- vector(length=n_obs_psa)
for(j in 1:n_obs_psa){
	psa.data$psa.age[j] <- psa.data$psa.time[j] + pt.data$age.dx[pt.data$id==psa.data$id[j]]}


psa.data$std.age <- (psa.data$psa.age-mean(psa.data$psa.age))/sd(psa.data$psa.age)
psa.data$age.basis <- get.ns.basis(obs.data=psa.data$std.age, knots=quantile(psa.data$std.age, p=c(0.025,0.5, 0.975)))

pt.data$std.vol <- rnorm(n,0,1)
pt.ordered <- pt.data[ordered,]
psa.data$std.vol <- vector(length=n_obs_psa)
for(j in 1:n_obs_psa){
	psa.data$std.vol[j] <- pt.data$std.vol[pt.data$id==psa.data$id[j]]}

b.vec <- matrix(nrow=n, ncol=3)
for(i in 1:n){
	b.vec[i,] <- mvrnorm(n=1, mu=mu.mat[,(pt.data$eta.true[i]+1)], Sigma=Sigma)}


apply(b.vec[pt.data$eta.true==0,],2,summary)
apply(b.vec[pt.data$eta.true==1,],2,summary)
write.csv(b.vec,"b-vec-true.csv")

psa.data$log.psa <- vector(length=n_obs_psa)
for(j in 1:n_obs_psa){
	lin.pred <- NULL
	lin.pred <- sum(b.vec[pt.data$id==psa.data$id[j],] * c(1, psa.data$std.age[j], psa.data$age.basis[j])) + beta[1]*psa.data$std.vol[j]
	if(pt.data$eta.true[pt.data$id==psa.data$id[j]]==1){
		lin.pred <- lin.pred + beta[2]*psa.data$std.vol[j]}
	psa.data$log.psa[j] <- rnorm(1,mean=lin.pred, sd=sigma_res)}


psa.data$subj<-vector(length=n_obs_psa)
for(j in 1:n_obs_psa){
	psa.data$subj[j] <- c(1:n)[ids.ordered=psa.data$id[j]]}


write.csv(psa.data,"psa-data-sim.csv")
write.csv(pt.data,"pt-data-sim.csv")
write.csv(pt.ordered,"pt-ordered-sim.csv")


