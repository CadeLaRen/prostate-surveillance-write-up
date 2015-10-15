library(dplyr)
library(ggplot2)

#IOP (Informative observation process)
IOP_BX <- TRUE #Informative observation process for biopsy
IOP_SURG <- TRUE #Informative observation process for surgery
leave_one_out <- FALSE

IOPs<-paste0(
	c('N')[!IOP_BX],'IOP_BX-',
	c('N')[!IOP_SURG],'IOP_SURG')
batch_path<-paste0(
	'batches/',
	IOPs,
	'/leave_one_out_',leave_one_out,
	'/batch-1/')
posterior_path<-paste0(batch_path,
	'concatenated_posterior.rds')

###############
# Load Data


psa_data_full<-read.csv("simulation-data/psa-data-sim.csv")
pt_data_full<-read.csv("simulation-data/pt-data-sim.csv") #true eta is now contained here, in patient data (pt data)
bx_data_full <- read.csv("simulation-data/bx-data-sim.csv")



### Output from leave-one-Out JAGS model  (output-out -> oo)
	# For now just use the output from the full model as an approximate replacement
### Then collect posterior estimates from full JAGS model
# of <- readRDS('batches/inf-obs/1/2015-07-24_JAGS_full_sample_summarized_12.rds')
oo <- readRDS(posterior_path)
of <- readRDS(posterior_path)



missing_etas <- which(is.na(pt_data_full$obs_eta)) #=1 if observed and aggressive, 0 if observed and not aggressive, or NA if not observed
eta_true_or_jags<-pt_data_full$obs_eta
eta_true_or_jags[missing_etas]<-colMeans(of$eta_hat_means) #indexed by subject

nsim <- dim(of1<-of$eta_hat_means)[1]
eta_jags1<-colMeans(of$eta_hat_means[1:nsim <  nsim/2,])
eta_jags2<-colMeans(of$eta_hat_means[1:nsim >= nsim/2,])
quantile(abs(eta_jags1-eta_jags2),prob=c(.5,.9,.99,1))
# From two MCMC samples the concordance is high. About 10 times higher than what we get with IS.

load(paste0(batch_path,'/2015-10-02_seed_101_IOP_BX-TRUE_IOP_SURG-TRUE_P-5e+06_online_fit_small_results.RData'))

if(FALSE){

files <- dir(batch_path)
fitfiles <- files[grep('online',files)]

small$draw_type <- 'small'
dynamic$draw_type <- 'dynamic'

bigConcatenated <- data.frame()
bigfiles <- files[grep('big',files)]
for(i in 1:length(bigfiles)){

	big_i <- readRDS(paste(batch_path,bigfiles[i],sep='/'))

}

}

	

big<-readRDS(paste0(batch_path,'2015-10-03_seed_101_IOP_BX-TRUE_IOP_SURG-TRUE_P-5e+06_online_fit_big_combined.rds'))
small$draw_type<-'small'
dynamic$draw_type<-'dynamic'
big$draw_type<-'big'


quantile(dynamic$time_IS,prob=c(.5,.9,.95,.99,1))

ofits<-rbind(small,dynamic,big) #out of sample fits

ofits$etas_jags<- eta_true_or_jags[ofits$subj]
ofits<-mutate(ofits, squared_errors_IS = (etas_jags-etas_IS)^2)
ofits<-mutate(ofits, errors_zhenke_abs = sqrt((etas_jags-etas_zhenke)^2) )

ofits<-mutate(ofits, errors_IS_abs = sqrt((etas_jags-etas_IS)^2) )
ofits<-mutate(ofits, errors_zhenke_abs = sqrt((etas_jags-etas_zhenke)^2) )


ggplot(ofits[ofits$draw_type=='small',]) +
	geom_point(aes(x=effective_ss, y=errors_IS_abs),cex=1.5) +
	scale_x_sqrt()+
	theme(text=element_text(size=18)) +
	labs(title='Absolute Error v. Calculation Time',x='Effective sample size\n(square root spacing)',y='Absolute difference', color='Candidate\ngenerating\nscheme')

ggplot(ofits) + geom_point(aes(x=time_IS, y=errors_IS_abs,color=factor(draw_type,levels=c('big','dynamic','small'))),cex=1.5) +
	scale_y_log10(breaks=c(.1,.05,.01,.005,.001)) +
	scale_x_log10(breaks=c(10^(1:7)))+
	theme(text=element_text(size=18)) +
	labs(title='Absolute Error v. Calculation Time',x='Calculation time for IS\n(seconds, log spacing)',y='Absolute difference (log spacing)', color='Candidate\ngenerating\nscheme')


ggplot(ofits) + geom_density(aes(x=time_IS, fill=draw_type),alpha=.6)+
	# scale_y_log10(breaks=c(.1,.05,.01,.005,.001)) +
	scale_x_log10(breaks=c(10^(1:7)))+
	theme(text=element_text(size=18)) +
	labs(title='Calculation Time',x='Calculation time for IS\n(seconds, log spacing)', color='Candidate\ngenerating\nscheme')

png(paste0('plots/',Sys.Date(),'_error_by_method.png'),height=400,width=460,pointsize=16)
ggplot(ofits) + geom_violin(aes(x=factor(draw_type,levels=c('small','dynamic','big')), y=errors_IS_abs),cex=1.5)+
	labs(title='Absolute Error v. Sampling Scheme',x='Candidate Generating Scheme',y='Absolute Difference')+
	theme(text=element_text(size=18))
dev.off()

# png(paste0('plots/',Sys.Date(),'_ESS_error_log.png'),height=400,width=460,pointsize=16)
ggplot(ofits) + geom_point(aes(x=effective_ss, y=errors_IS_abs,color=draw_type),cex=0.8) +
	scale_y_log10(breaks=c(.1,.05,.01,.005,.001)) +
	scale_x_log10(breaks=c(10^c(1:7)))+
	theme(text=element_text(size=18)) +
	labs(title='Effective Sample Size v. Absolute Error',x='Effective sample size for IS (log spacing)',y='Absolute difference (log spacing)', color='Candidate\ngenerating\nscheme') +
	geom_vline(xintercept=1000,lty=2)+
	theme(axis.text.x=element_text(angle=35, hjust = 1))
dev.off()

ggplot(ofits[ofits$draw_type!='big',]) + geom_line(aes(x=effective_ss, y=errors_IS_abs, group=subj)) +
	# scale_y_log10(breaks=c(.1,.05,.01,.005,.001)) +
	scale_x_log10(breaks=c(10^(1:7)))+
	labs(title='Effective Sample Size v. Absolute Error',x='Effective sample size for IS (log spacing)',y='Absolute difference (log spacing)', color='Candidate\ngenerating\nscheme') +
	geom_vline(xintercept=1000,lty=2)

group_by(ofits,draw_type) %>%
	summarize(mean(errors_IS_abs))


ggplot(ofits) + geom_point(aes(x=effective_ss, y=squared_errors_IS,color=draw_type),cex=0.9) +
	geom_vline(xintercept=1000,lty=2)+
	theme(text=element_text(size=18)) +
	labs(title='Effective Sample Size v. Absolute Error',x='Effective sample size for IS\n(sqrt spacing)',y='Absolute difference', color='Candidate\ngenerating\nscheme')+
	theme(axis.text.x=element_text(angle=35, hjust = 1))


# png(paste0('plots/',Sys.Date(),'_ESS_error_sqrt.png'),height=440,width=450,pointsize=16)
ggplot(ofits) + geom_point(aes(x=effective_ss, y=errors_IS_abs,color=draw_type),cex=0.9) +
	scale_x_sqrt()+
	geom_vline(xintercept=1000,lty=2)+
	theme(text=element_text(size=18)) +
	labs(title='Effective Sample Size v. Absolute Error',x='Effective sample size for IS\n(sqrt spacing)',y='Absolute difference', color='Candidate\ngenerating\nscheme')+
	theme(axis.text.x=element_text(angle=35, hjust = 1))
dev.off()

ggplot(ofits[ofits$draw_type=='small',]) +
	geom_point(aes(x=sqrt(effective_ss), y=errors_IS_abs),cex=1.5) +
	# scale_y_log10(breaks=c(.1,.05,.01,.005,.001)) +
	# scale_x_log10(breaks=c(10^(1:7)))+
	theme(text=element_text(size=18)) +
	labs(title='Effective Sample Size v. Absolute Error',x='sqrt(effective sample size for IS))',y='Absolute difference')

# png(paste0('plots/',Sys.Date(),'_agreement_IS_MCMC.png'),height=400,width=450,pointsize=16)
ggplot(ofits)+geom_point(aes(x=etas_IS,y=etas_jags,color=draw_type),cex=1.2) +
	theme(text=element_text(size=18)) +
	labs(title='Agreement of risk estimates from\nIS and MCMC',x='Risk estimates from IS',y='Risk estimates from MCMC', color='Candidate\ngenerating\nscheme')
dev.off()



