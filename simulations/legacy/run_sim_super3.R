#########################################################################################
#Run Simulation
##########################################################################################
rm(list = ls())
gc()
path <- "~/Spike_Slab/Final/"
path2 <- "~/Spike_Slab/Final/Superinduction_Results/"
source(paste0(path,"EP_Functions2.R"))
source(paste0(path,"support_functions2.R"))
source("~/Spike_Slab/Final/FASS_Sampler4.R")
source("~/Spike_Slab/Final/VB_functions.R")
source("~/Spike_Slab/Final/Laplace_Approx.R")



p.seq <- c(500)
n.seq <- c(5000,10000,25000,50000,75000,100000,150000,200000,250000,300000,400000,500000,750000,10^6)
#n.seq <- c(5000,10000,25000,50000,100000,250000,500000,10^6)
p.true <- c(10)
L.t <- length(p.true)
L.p <- length(p.seq)
L.n <- length(n.seq)
col.num <- L.t*L.p*L.n
times <- matrix(nrow=5,ncol=col.num) -> TPR
TNR <- matrix(nrow=5,ncol=col.num)
#for(k1 in 1:length(n.seq)){
#for(k2 in 1:length(p.seq)){
for(k3 in 1:length(p.true)){
  #set.seed(4747)
  n <- max(n.seq)
  p <- max(p.seq)
  true.num <- p.true[k3] 
  J = matrix(0,nrow=n,ncol=p)
  
  set.seed(4747)
  X <- rep(1,n)
  K <- matrix(0,nrow=n,ncol=true.num)
  means <- rep(0,true.num)
  sd <- rchisq(true.num,20)
  for(i in 1:true.num) K[,i] <- rnorm(n,means[i],sd[i])
  print("K Made")
  alpha <- rnorm(1,0,sqrt(10))
  beta <- rnorm(true.num,0,10)
  
  mu_delta <- rnorm(n,0,sqrt(10))
  y <- rnorm(n,X*alpha+K%*%beta[1:true.num]+mu_delta,1)
  epsilon <- 0.01
  
  
  for(i in 1:true.num) J[,i] <- K[,i]
  
  for(i in (true.num+1):p){
    J.temp <- (1-epsilon)*rnorm(n,0,sqrt(10)) + epsilon*mu_delta
    J[,i] <- J.temp
    if(i%%20==0) print(round(100*i/p,2)) }

rm(K)
 
J <- scale(J, center = FALSE, scale = TRUE)
scales <- attr(J, "scaled:scale")
saveRDS(scales, file=paste0(path2,'scales_VB',collapse='_'))  
  print("data made")
  
  

  for(k1 in 1:length(n.seq)){
    n <- n.seq[k1]
    for(k2 in 1:L.p){
      p <- p.seq[k2]
      start <- Sys.time()
      gamma.out_laplace <- FASS.sampler(y[1:n],J[1:n,1:p],10000, 1, 1,"Laplace",1, sqrt(n),full_posterior='Block-VB')
      times[1,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- difftime(Sys.time(),start,units="secs")
      TPR[1,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(gamma.out_laplace$gamma[which.max(gamma.out_laplace$odds)-1,1:true.num])
      TNR[1,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- 1- mean(gamma.out_laplace$gamma[which.max(gamma.out_laplace$odds)-1,(true.num+1):p])
      saveRDS(gamma.out_laplace, file=paste0(path2,'_results_SSL_VB_',paste0(c(n,p,true.num),collapse='_')))
      rm(gamma.out_laplace)
      print(Sys.time()-start)
      start <- Sys.time()
      
      gamma.out_gauss <- FASS.sampler(y[1:n],J[1:n,1:p],10000, 1, 1,"Gaussian",1, n,full_posterior='Block-VB')
      times[2,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- difftime(Sys.time(),start,units="secs")
      TPR[2,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(gamma.out_gauss$gamma[which.max(gamma.out_gauss$odds)-1,1:true.num])
      TNR[2,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- 1- mean(gamma.out_gauss$gamma[which.max(gamma.out_gauss$odds)-1,(true.num+1):p])
      print(Sys.time()-start)
      saveRDS(gamma.out_gauss, file=paste0(path2,'_results_SSG_VB_',paste0(c(n,p,true.num),collapse='_')))  
rm(gamma.out_gauss)

    
      start <- Sys.time()
      gamma.out_student1 <- FASS.sampler(y[1:n],J[1:n,1:p],10000, 1, 1,"Student",1, n,2,full_posterior='Block-VB')
      times[3,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- difftime(Sys.time(),start,units="secs")
      TPR[3,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(gamma.out_student1$gamma[which.max(gamma.out_student1$odds)-1,1:true.num])
      TNR[3,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- 1- mean(gamma.out_student1$gamma[which.max(gamma.out_student1$odds)-1,(true.num+1):p])
      print(Sys.time()-start)
      saveRDS(gamma.out_student1, file=paste0(path2,'_results_SSS1_VB_',paste0(c(n,p,true.num),collapse='_')))
      rm(gamma.out_student1)
      
      start <- Sys.time()
      gamma.out_student2 <- FASS.sampler(y[1:n],J[1:n,1:p],10000, 1, 1,"Student",1, n,5,full_posterior='Block-VB')
      times[4,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- difftime(Sys.time(),start,units="secs")
      TPR[4,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(gamma.out_student2$gamma[which.max(gamma.out_student2$odds)-1,1:true.num])
      TNR[4,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- 1- mean(gamma.out_student2$gamma[which.max(gamma.out_student2$odds)-1,(true.num+1):p])
      print(Sys.time()-start)
      saveRDS(gamma.out_student2, file=paste0(path2,'_results_SSS2_VB_',paste0(c(n,p,true.num),collapse='_')))
      rm(gamma.out_student2)
      
      start <- Sys.time()
      gamma.out_student3 <- FASS.sampler(y[1:n],J[1:n,1:p],10000, 1, 1,"Student",1, n,10,full_posterior='Block-VB')
      times[5,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- difftime(Sys.time(),start,units="secs")
      TPR[5,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(gamma.out_student3$gamma[which.max(gamma.out_student3$odds)-1,1:true.num])
      TNR[5,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- 1- mean(gamma.out_student3$gamma[which.max(gamma.out_student3$odds)-1,(true.num+1):p])
      print(Sys.time()-start)
      saveRDS(gamma.out_student3, file=paste0(path2,'_results_SSS3_VB_',paste0(c(n,p,true.num),collapse='_')))
      rm(gamma.out_student3)
      
      print(c(n,p,true.num))
      print(times[,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)])
      print(TPR[,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)])
      print(TNR[,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)])
      
      saveRDS(times, file=paste0(path2,'_results_times_VB'))
      saveRDS(TPR, file=paste0(path2,'_results_TPR_VB'))
      saveRDS(TNR, file=paste0(path2,'_results_TNR_VB'))
      
    }}}
print(times)


