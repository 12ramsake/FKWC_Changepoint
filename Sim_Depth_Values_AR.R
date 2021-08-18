




#This script is for the AR simulation in the paper.
#It only computes the depth values for each of the 200 runs . 
#You must then run the AMOC or PELT alorithm on the ranks of the depth values after





#you may not need all of these packages
library(sn)
library(fda)
library(roahd)
library(mvtnorm)
library(fdasrvf)
library(doParallel)
library(MFHD)
library(mrfDepth)
library(doRNG)
#For this script you need the function gen_AR, 
#which I built using code from Dr. Martin Wendler, since it is not my code I have not shared it here
#You should see the paper Bootstrapping covariance operators of functional time series by Sharipov et. al 
#and use that as a template to build a gen_AR function 


FMp=function(data,derivs){
  
  dp=sapply(1:100,function(x){ddalpha::depth.halfspace(cbind( data[,x],derivs[,x]),cbind( data[,x],derivs[,x]),num.directions =100)   })
  return(rowMeans(dp))
}


RPDd=function(data,derivs){
  gen_direction=replicate(20,rnorm(100))
  gen_direction=apply(gen_direction,2,function(x){x/sqrt(sum(x^2))})
  
  #data is n by m
  # dim(data)
  bd=apply(gen_direction,2,function(u){
    bdd=cbind(data%*%u,derivs%*%u)
    ddalpha::depth.simplicial(bdd,bdd)})
  return(rowMeans(bd))
}


RPD=function(data){
  gen_direction=replicate(20,rnorm(100))
  gen_direction=apply(gen_direction,2,function(x){x/sqrt(sum(x^2))})
  
  #data is n by m
  # dim(data)
  bd=apply(gen_direction,2,function(u){
    bdd=data%*%u
    rnk=rank(bdd)
    rnk*(length(bdd)-rnk)})
  
  return(rowMeans(bd))
}

#simulate data and the resulting depths using the derivative
simulate_mv_depths<-function(N1,N2,D1,D2){
  
  #columns are data, so transpose
  fdat=t(gen_AR(N1+N2,D1,D2))
  dat=list(mdata=fdat,argvals=seq(0,1,l=100))
  rownames(dat$mdata)=as.character(1:(2*N1))
  fdat<-fdata(dat$mdata,dat$argvals)
  
  derivs<-derivcurves(dat$mdata)
  fderivs<-fdata(derivs,dat$argvals)
  
  RPD_depth=RPD(dat$mdata)
  FM_depth=depth.FM(fdat)$dep
  LTR_depth=c(norm.fdata(fdat))
  
  
  RPD_depth_d=RPDd(dat$mdata,derivs)
  FM_depth_d=FMp(dat$mdata,derivs)
  LTR_depth_d=c(norm.fdata(fdat))+c(norm.fdata(fderivs))
  
  
  depths=data.frame(FM_depth,RPD_depth,LTR_depth,FM_depth_d,RPD_depth_d,LTR_depth_d)
  
  return(depths)
  
}



runMVSim<-function(N1,N2,D1,D2,num_runs,fileName){
  
  
  no_cores<-detectCores()-1
  
  cl <- makeCluster(no_cores,type="FORK")   
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  
  depth_values=try({foreach(i=1:num_runs) %dopar% {simulate_mv_depths(N1,N2,D1,D2)}})
  
  errorsp=inherits(depth_values, "try-error")
  if(errorsp){
    print("there was an error!")
    
  }
  
  stopCluster(cl)
  registerDoSEQ()
  closeAllConnections()
  
  #Your simulation folder
  dirr=""
  save(depth_values,file=paste(dirr,fileName,sep=""))
  # }
}






num_runs=200
grid=seq( 0, 1, length.out = 100 )


N1=N2=50

d_mat=rbind(c(0,0),
            c(0.4,0),
            c(0.8,0),
            c(0,0.4),
            c(0,0.8),
            c(0.4,0.4))


for(i in 1:nrow(d_mat)){
  d1=d_mat[i,1]
  d2=d_mat[i,2]
  fileName=paste0("AR_Shar_both_depths_N_100_d1_",d1,"_d2_",d2,"_num_runs_200.Rda",sep="")
  runMVSim(N1,N2,d1,d2,num_runs,fileName)
  
}







