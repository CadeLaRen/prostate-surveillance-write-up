



#All values passed to `args` are logical
args<-as.logical(commandArgs(TRUE))
	#The base case is:
	# args = c(TRUE, TRUE, FALSE)

IOP_BX <- args[1]
IOP_SURG <- args[2]
leave_one_out <- args[3]

taskID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if(is.na(taskID)) taskID <- 0

#import environment variable, used for running multiple chains in parallel
if( leave_one_out) {SEED <- 0		; star <- taskID; 	crop <- TRUE}
if(!leave_one_out) {SEED <- taskID 	; star <- 0; 		crop <- FALSE} #here, parallel is just used to extend the number of MCMC draws.




(save_path <- paste0(
	'batches/',
	c('N')[!IOP_BX],'IOP_BX-',
	c('N')[!IOP_SURG],'IOP_SURG',
	'/leave_one_out_',leave_one_out,'/batch-1/'))

if(!dir.exists(save_path)) dir.create(save_path,recursive=TRUE)

# MCMC settings
#ni <- 100; nb <- 20; nt <- 5; nc <- 1
# ni <- 1000; nb <- 20; nt <- 5; nc <- 1
ni <- 50000; nb <- 25000; nt <- 20; nc <- 1
# ni <-  250000; nb <- 25000; nt <- 10; nc <- 1 ##Don't use this setting!! It results in vectors too large for memory and often crashes. Uses around 4G vmem, and about 30 hours. 
	# Problem could be because we're storing all those subject specific variables
	#JAGS call often spits out: 
		#terminate called after throwing an instance of 'std::bad_alloc'
  		# what():  std::bad_alloc
	

(len.sim <- round((ni - nb)/nt))

cat('',file=paste0(save_path,Sys.Date(),'_started_SEED-',SEED,'_star-',star,'_crop-',crop,'_len.sim-',len.sim,'.txt'))

source('call-jags-master.R') #returns function `do_one`
outJAGS <- do_one(seed=SEED) #final output of this script


#Reduce filesize of results (make new object, outJAGS2)
#Parameters to ignore
params2ignore<-c('b_vec','rho_int_slope','eta_hat','p_rc')
if(IOP_BX) params2ignore <- unique(c(params2ignore, 'p_bx'))
if(IOP_SURG) params2ignore <- unique(c(params2ignore, 'p_surg'))


outJAGS2 <- outJAGS[!names(outJAGS) %in% c('sims.matrix','sims.list','sims.array')]

#Now define the sims.list for outJAGS2
outJAGS2_sl <- outJAGS$sims.list[!names(outJAGS$sims.list) %in% params2ignore]
outJAGS2_sl$eta_hat_means <- apply(outJAGS$sims.list$eta_hat,2,mean)
if(crop){ #if `crop`, save subj-specific posterior for subj star. otherwise don't.
	outJAGS2_sl$b_vec.star <- outJAGS$sims.list$b_vec[,star,]
	outJAGS2_sl$eta_hat.star <- NULL
	if(star > 214)
		outJAGS2_sl$eta_hat.star <- outJAGS$sims.list$eta_hat[,star - 214]
}

outJAGS2$sims.list <- outJAGS2_sl


corner <- function(x, n=6 ){
	d <- dim(x)
	ld<-length(d)

	if(ld==0) return(x[1:min(n,length(x))])

	#If x is an array:
	text_n <- paste0('1:',pmin(n,dim(x)))
	text_ind <- paste0('x[', paste(text_n,collapse=','), ']')
	eval(parse(text=text_ind))
}
sizeMb<-function(x) {object.size(x)/1000000}

lapply(outJAGS$sims.list,corner)
lapply(outJAGS2$sims.list,sizeMb)
sizeMb(outJAGS2)

saveRDS(outJAGS2,file=paste0(save_path,Sys.Date(),'_posterior_IOP_BX-',IOP_BX,'_IOP_SURG-',IOP_SURG,'_SEED-',SEED,'_star-',star,'_crop-',crop,'_nsim-',len.sim,'.rds'))
saveRDS(lapply(outJAGS$sims.list,corner),file=paste0(save_path,Sys.Date(),'_full_posterior_corner_IOP_BX-',IOP_BX,'_IOP_SURG-',IOP_SURG,'_SEED-',SEED,'_star-',star,'_crop-',crop,'_nsim-',len.sim,'.rds')) #note, this is *not* outJAGS2.
str(outJAGS2$sims.list)
summary(outJAGS2$sims.list$mu_slope)

