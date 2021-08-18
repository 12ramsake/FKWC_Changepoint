

# This script is for the main AMOC simulation in the paper.
#It only computes the depth values for each of the 200 runs . 
#You must then run the AMOC or PELT algorithm on the ranks of the depth values after


library(sn)
library(fda)
library(roahd)
library(mvtnorm)
library(fdasrvf)
library(doParallel)
library(MFHD)
library(mrfDepth)
library(doRNG)
####kernels and simulating data

#beta will be the scale parameter for any kernel
#different kernels
k_linear<-function(t1,t2,beta,c,inte){inte+beta*(t1-c)*(t2-c)}
k_exp<-function(t1,t2,beta,alpha,dummy){beta*exp(-(t1-t2)^2/(2*alpha^2))}
# k_exp2<-function(t1,t2,beta,alpha){alpha*exp(-beta*abs(t1-t2))}
k_periodic<-function(t1,t2,beta,alpha,p){beta*exp(-2/alpha^2*(sin(pi*(t1-t2)/p)^2))}


##get covariance over a grid, makes pd if not bc rounding
getCov<-function(grid,k){
  return(corpcor::make.positive.definite(outer(grid,grid,k)))
}

###generate a set of functional data with specific kernel, 
##and N is sample size

gfd<-function(N1=10,N2=10,
              res=1e2, 
              Cov1,
              Cov2,
              dist="t",
              del=0.9){
  
  
  grid = seq( 0, 1, length.out = res )
  
  
  if(dist=="N"){
    Data1 = generate_gauss_fdata( N1,  rep(0,res),Cov = Cov1 )
    Data2 = generate_gauss_fdata( N2, rep(0,res), Cov = Cov2 )
  }
  else if(dist=="t"){
    Data1 = rmvt( N1,  sigma =  Cov1/3,df=3 )
    Data2 = rmvt( N2,  sigma =  Cov2/3,df=3 )
  }
  else{
    cp1=list(mean=rep(0,nrow(Cov1)), var.cov=Cov1, gamma1=rep(del,nrow(Cov1))/nrow(Cov1))
    cp2=list(mean=rep(0,nrow(Cov2)), var.cov=Cov2, gamma1=rep(del,nrow(Cov2))/nrow(Cov2))
    dp1=cp2dp(cp1, "SN")
    dp2=cp2dp(cp2, "SN")
    
    Data1 = rmsn(N1, dp=dp1)
    Data2 = rmsn(N1, dp=dp2)
  }
  
  return(list(argvals =  grid, mdata=rbind(Data1,Data2)))
  
}




runMVSim<-function(N1,N2,c1,c2,grid,num_runs,fileName,
                   warpF=F,warp_v=1,warp_all=F,outlier_type=0,
                   dist="t",del=0.9,windows=F){
  
  
  no_cores<-detectCores()-1
  
  cl <- makeCluster(no_cores,type="FORK")   
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  
  depth_values=try({foreach(i=1:num_runs) %dopar% {simulate_mv_depths(N1,N2,c1,c2,res=100,del=del,dist=dist,grid=grid,warpF=warpF,warp_sigma=warp_v,outlier_type=outlier_type,warp_all=warp_all)}})
  
  errorsp=inherits(depth_values, "try-error")
  if(errorsp){
    print("there was an error!")
    
  }
  
  stopCluster(cl)
  registerDoSEQ()
  closeAllConnections()
  
  #your sim directory
  dirr=""
  save(depth_values,file=paste(dirr,fileName,sep=""))
  # }
}




#Source the initial script J2
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
simulate_mv_depths<-function(N1,N2,c1,c2,res=100,
                             grid=seq( 0, 1, length.out = res ),
                             dist="t",
                             del=0.9,
                             warpF=F,
                             warp_all=F,
                             warp_sigma=1,
                             outlier_type=0){
  
  dat=gfd(N1,N2,Cov1=c1,Cov2=c2,res=res,dist,del)
  rownames(dat$mdata)=as.character(1:(2*N1))
  fdat<-fdata(dat$mdata,dat$argvals)
  
  derivs<-derivcurves(dat$mdata)
  fderivs<-fdata(derivs,dat$argvals)
  
  RPD_depth=RPD(dat$mdata)
  FM_depth=depth.FM(fdat)$dep
  LTR_depth=c(norm.fdata(fdat))
  
  
  RPD_depth_d=RPDd(dat$mdata,derivs)
  FM_depth_d=FMp(dat$mdata,derivs)
  
  #LTR
  LTR_depth_d=c(norm.fdata(fdat))+c(norm.fdata(fderivs))
  
  
  depths=data.frame(FM_depth,RPD_depth,LTR_depth,FM_depth_d,RPD_depth_d,LTR_depth_d)
  
  return(depths)
  
}


num_runs=200
grid=seq( 0, 1, length.out = 100 )



for(N1 in c(50,100,250)){
  
  N2=N1
  print("sim 1")
  #changes in beta
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  tmp=seq(0.5,1,length.out = 15)
  for(dist in c('N','t','sn')){
    for(betai in 1:15){
      beta=tmp[betai]
      fileName=paste0("FKWC_Results_cp/dist_",dist,"_both_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,beta,.05)}
      c2=getCov(grid,k2)
      runMVSim(N1,N2,c1,c2,grid,num_runs,fileName,dist="t")
      count=count+1
      print(count/(15*3))
    }
  }
  
  
  
  print("sim 2")
  #changes in alpha
  count=0
  k1=function(t1,t2){k_exp(t1,t2,0.5,.05)}
  c1=getCov(grid,k1)
  tmp=seq(0.05,0.1,length.out=15)
  for(dist in c('N','t','sn')){
    for(alphai in 1:15){
      alpha=tmp[alphai]
      fileName=paste0("FKWC_Results_cp/dist_",dist,"_both_depths_k_exp_alpha1_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
      c2=getCov(grid,k2)
      runMVSim(N1,N2,c1,c2,grid,num_runs,fileName,dist="t")
      count=count+1
      print(count/(15*3))
    }
  }
  
  
}











