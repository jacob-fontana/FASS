#########################################################################################
#Run Simulation
##########################################################################################
rm(list = ls())
gc()
path <- "~/Spike_Slab/Final/"
path2 <- "~/Spike_Slab/Final/Superinduction_Results/"
source("~/Spike_Slab/Final/support_functions2.R")
source("~/Spike_Slab/Final/Gibbs_samplers_cpp.R")



p.seq <- c(500)
#n.seq <- c(5000,10000,25000,50000,75000,100000,150000,200000,250000,300000,400000,500000,750000,10^6)
n.seq <- c(5000,10000,25000,50000,100000,250000,500000,10^6)
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
saveRDS(scales, file=paste0(path2,'scales_Skinny',collapse='_'))    

  
  print("data made")
  
  

  for(k1 in 1:length(n.seq)){
    n <- n.seq[k1]
    for(k2 in 1:L.p){
      p <- p.seq[k2]
      start <- Sys.time()

      gamma.out_laplace <- skinnyGibbs(J[1:n,1:p], y[1:n], 1/n, 1, 1,1, nburn=2000, niter=8000,family="Laplace")
      times[1,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- difftime(Sys.time(),start,units="secs")
      TPR[1,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(colMeans(gamma.out_laplace$outgamma[,1:true.num])>0.5)
      TNR[1,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(colMeans(gamma.out_laplace$outgamma[,(true.num+1):p])<0.5)
      saveRDS(gamma.out_laplace, file=paste0(path2,'_results_SSL_Skinny_',paste0(c(n,p,true.num),collapse='_')))
      print(Sys.time()-start)
      start <- Sys.time()
      
      gamma.out_gauss <- skinnyGibbs(J[1:n,1:p], y[1:n], 1/n, 1, 1,1, nburn=2000, niter=8000,family="Gaussian")
      times[2,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- difftime(Sys.time(),start,units="secs")
      TPR[2,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(colMeans(gamma.out_gauss$outgamma[,1:true.num])>0.5)
      TNR[2,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(colMeans(gamma.out_gauss$outgamma[,(true.num+1):p])<0.5)
      print(Sys.time()-start)
      saveRDS(gamma.out_gauss, file=paste0(path2,'_results_SSG_Skinny_',paste0(c(n,p,true.num),collapse='_')))  
    
      start <- Sys.time()
      gamma.out_student1 <- skinnyGibbs(J[1:n,1:p], y[1:n], 1/n, 1, 1,1, nburn=2000, niter=8000,family="Student",nu=2)
      times[3,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- difftime(Sys.time(),start,units="secs")
      TPR[3,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(colMeans(gamma.out_student1$outgamma[,1:true.num])>0.5)
      TNR[3,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(colMeans(gamma.out_student1$outgamma[,(true.num+1):p])<0.5)
      print(Sys.time()-start)
      saveRDS(gamma.out_student1, file=paste0(path2,'_results_SSS1_Skinny_',paste0(c(n,p,true.num),collapse='_')))
      
      start <- Sys.time()
      gamma.out_student2 <- skinnyGibbs(J[1:n,1:p], y[1:n], 1/n, 1, 1,1, nburn=2000, niter=8000,family="Student",nu=5)
      times[4,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- difftime(Sys.time(),start,units="secs")
      TPR[4,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(colMeans(gamma.out_student2$outgamma[,1:true.num])>0.5)
      TNR[4,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(colMeans(gamma.out_student2$outgamma[,(true.num+1):p])<0.5)
      print(Sys.time()-start)
      saveRDS(gamma.out_student2, file=paste0(path2,'_results_SSS2_Skinny_',paste0(c(n,p,true.num),collapse='_')))
      
      start <- Sys.time()
      gamma.out_student3 <- skinnyGibbs(J[1:n,1:p], y[1:n], 1/n, 1, 1,1, nburn=2000, niter=8000,family="Student",nu=10)
      times[5,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- difftime(Sys.time(),start,units="secs")
      TPR[5,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(colMeans(gamma.out_student3$outgamma[,1:true.num])>0.5)
      TNR[5,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)] <- mean(colMeans(gamma.out_student3$outgamma[,(true.num+1):p])<0.5)
      print(Sys.time()-start)
      saveRDS(gamma.out_student3, file=paste0(path2,'_results_SSS3_Skinny_',paste0(c(n,p,true.num),collapse='_')))
      rm(gamma.out_student3)
      
      print(c(n,p,true.num))
      print(times[,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)])
      print(TPR[,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)])
      print(TNR[,((k1-1)*L.t*L.p+(k2-1)*L.t+k3)])
      
      saveRDS(times, file=paste0(path2,'_results_times_Skinny'))
      saveRDS(TPR, file=paste0(path2,'_results_TPR_Skinny'))
      saveRDS(TNR, file=paste0(path2,'_results_TNR_Skinny'))
      
    }}}
print(times)