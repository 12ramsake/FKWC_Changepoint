
# This script is for the main  5 change-point simulation in the paper.
#It only computes the depth values for each of the 200 runs . 
#You must then run the  PELT algorithm on the ranks of the depth values after
#You can use the FKWC Methods script to do this 



library(sn)
library(MFHD)
library(mrfDepth)
library(fda)
library(fda.usc)
library(roahd)
library(mvtnorm)
# library(fdasrvf)
library(doParallel)
library(doRNG)




# at least 10% apart
compute_location_change_points=function(N){
  min_length=N*.1
  candidate=sort(runif(5,.1,.9))
  while(sum(diff(candidate)<.1)>0)
    candidate=sort(runif(5,.1,.9))
  changes=floor(candidate*N)
  N1=changes[1]
  N2=changes[2]-changes[1]
  N3=changes[3]-changes[2]
  N4=changes[4]-changes[3]
  N5=changes[5]-changes[4]
  N6=N-changes[5]
  return(c(N1,N2,N3,N4,N5,N6))
}



set.seed(440)
# changes_100=replicate(num_runs, compute_location_change_points(100))
changes_200=replicate(num_runs, compute_location_change_points(200))
changes_500=replicate(num_runs, compute_location_change_points(500))
changes_1000=replicate(num_runs, compute_location_change_points(1000))
changes_2500=replicate(num_runs, compute_location_change_points(2500))
# changes=changes/1000
#reset
set.seed(NULL)













####kernels and simulating data

k_exp<-function(t1,t2,beta,alpha,dummy){beta*exp(-(t1-t2)^2/(2*alpha^2))}


##get covariance over a grid, makes pd if not bc rounding
getCov<-function(grid,k){
  return(corpcor::make.positive.definite(outer(grid,grid,k)))
}



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







###generate a set of functional data with specific kernel, 
##and N is sample size
#generate data from spec. covariance kernel
gfd_k_5<-function(sample_sizes,
                  Cov1,
                  Cov2,
                  Cov3,
                  Cov4,
                  Cov5,
                  Cov6,
                  dist="t",
                  del=0.9,
                  res=1e2){
  
  
  grid = seq( 0, 1, length.out = res )
  # sample_sizes=compute_location_change_points(N)
  
  # for(j in 1:4){}
  if(dist=="N"){
    Data1 = generate_gauss_fdata( sample_sizes[1],  rep(0,res),Cov = Cov1 )
    Data2 = generate_gauss_fdata( sample_sizes[2], rep(0,res), Cov = Cov2 )
    Data3 = generate_gauss_fdata( sample_sizes[3], rep(0,res), Cov = Cov3 )
    Data4 = generate_gauss_fdata( sample_sizes[4], rep(0,res), Cov = Cov4 )
    Data5 = generate_gauss_fdata( sample_sizes[5], rep(0,res), Cov = Cov5 )
    Data6 = generate_gauss_fdata( sample_sizes[6], rep(0,res), Cov = Cov6 )
  }
  else if(dist=="t"){
    Data1 = rmvt( sample_sizes[1],  sigma =  Cov1/3, df=3 )
    Data2 = rmvt( sample_sizes[2],  sigma =  Cov2/3, df=3 )
    Data3 = rmvt( sample_sizes[3],  sigma =  Cov3/3 , df=3 )
    Data4 = rmvt( sample_sizes[4],  sigma =  Cov4/3, df=3 )
    Data5 = rmvt( sample_sizes[5],  sigma =  Cov5/3, df=3 )
    Data6 = rmvt( sample_sizes[6],  sigma =  Cov6/3, df=3 )
  }
  else{
    # Data1 = rdmsn(N1, res, rep(0,res),Cov1, del=rep(del,res))
    # Data2 = rdmsn(N2, res, rep(0,res),Cov2, del=rep(del,res))
    # Data3 = rdmsn(N3, res, rep(0,res),Cov3, del=rep(del,res))
    
    
    cp1=list(mean=rep(0,nrow(Cov1)), var.cov=Cov1, gamma1=rep(del,nrow(Cov1))/nrow(Cov1))
    cp2=list(mean=rep(0,nrow(Cov2)), var.cov=Cov2, gamma1=rep(del,nrow(Cov2))/nrow(Cov2))
    cp3=list(mean=rep(0,nrow(Cov3)), var.cov=Cov3, gamma1=rep(del,nrow(Cov3))/nrow(Cov3))
    cp4=list(mean=rep(0,nrow(Cov4)), var.cov=Cov4, gamma1=rep(del,nrow(Cov4))/nrow(Cov4))
    cp5=list(mean=rep(0,nrow(Cov5)), var.cov=Cov5, gamma1=rep(del,nrow(Cov5))/nrow(Cov5))
    cp6=list(mean=rep(0,nrow(Cov6)), var.cov=Cov6, gamma1=rep(del,nrow(Cov6))/nrow(Cov6))
    
    dp1=cp2dp(cp1, "SN")
    dp2=cp2dp(cp2, "SN")
    dp3=cp2dp(cp3, "SN")
    dp4=cp2dp(cp4, "SN")
    dp5=cp2dp(cp5, "SN")
    dp6=cp2dp(cp6, "SN")
    
    Data1 = rmsn(sample_sizes[1], dp=dp1)
    Data2 = rmsn(sample_sizes[2], dp=dp2)
    Data3 = rmsn(sample_sizes[3], dp=dp3)
    Data4 = rmsn(sample_sizes[4], dp=dp4)
    Data5 = rmsn(sample_sizes[5], dp=dp5)
    Data6 = rmsn(sample_sizes[6], dp=dp6)
    
  }
  
  
  
  
  
  
  
  
  return(list(argvals =  grid, mdata=rbind(Data1,Data2,Data3,Data4,Data5,Data6)))
  
}




#simulate epidemic model and then compute depth functions
simulate_mv_depths<-function(sample_sizes,
                             Cov1,
                             Cov2,
                             Cov3,
                             Cov4,
                             Cov5,
                             Cov6,
                             res=100,
                             dist="t",
                             del=0.9){
  
  dat=gfd_k_5(sample_sizes,Cov1,Cov2,Cov3,Cov4,Cov5,Cov6,res=res,dist=dist,del=del)
  rownames(dat$mdata)=as.character(1:sum(sample_sizes))
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




#simulate data and the resulting depths using the derivative
runMVSim<-function(N=50,
                   Cov1,
                   Cov2,
                   Cov3,
                   Cov4,
                   Cov5,
                   Cov6,
                   grid,num_runs,
                   fileName,
                   dist="t",
                   del=0.9){
  
  
  
  no_cores<-75
  
  cl <- makeCluster(no_cores,type="FORK")   
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  
  if(N==200)
    depth_values=try({foreach(i=1:num_runs) %dopar% {simulate_mv_depths(changes_200[,i],Cov1,Cov2,Cov3,Cov4,Cov5,Cov6,res=1e2,del=del,dist=dist)}})
  if(N==500)
    depth_values=try({foreach(i=1:num_runs) %dopar% {simulate_mv_depths(changes_500[,i],Cov1,Cov2,Cov3,Cov4,Cov5,Cov6,res=1e2,del=del,dist=dist)}})
  if(N==1000)
    depth_values=try({foreach(i=1:num_runs) %dopar% {simulate_mv_depths(changes_1000[,i],Cov1,Cov2,Cov3,Cov4,Cov5,Cov6,res=1e2,del=del,dist=dist)}})
  if(N==2500)
    depth_values=try({foreach(i=1:num_runs) %dopar% {simulate_mv_depths(changes_2500[,i],Cov1,Cov2,Cov3,Cov4,Cov5,Cov6,res=1e2,del=del,dist=dist)}})
  
  
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
}





##


num_runs=200
grid=seq( 0, 1, length.out = 100 )





for(N in c(200,500,1000,2500)){
  # N2=N1
  # N3=N1
  print("sim 1")
  #changes in beta
  count=0
  #baseline
  k1=function(t1,t2){k_exp(t1,t2,0.5,0.05)}
  c1=getCov(grid,k1)
  k2=function(t1,t2){k_exp(t1,t2,1,0.05)}
  c2=getCov(grid,k2)
  k3=function(t1,t2){k_exp(t1,t2,.5,0.05)}
  c3=getCov(grid,k3)
  k4=function(t1,t2){k_exp(t1,t2,1,0.05)}
  c4=getCov(grid,k4)
  k5=function(t1,t2){k_exp(t1,t2,.5,0.05)}
  c5=getCov(grid,k5)
  k6=function(t1,t2){k_exp(t1,t2,1,0.05)}
  c6=getCov(grid,k6)
  
  # tmp=seq(0.5,1,length.out = 15)
  for(dist in c("N","t","sn")){
    fileName=paste0("FKWC_Results_cp/dist_",dist,"_mfn_depths_k_exp_Beta_alt_mult_N_",N,"_num_runs_200.Rda",sep="")
    runMVSim(N,c1,c2,c3,c4,c5,c6,grid,num_runs,fileName,dist=dist)
    count=count+1
    print(count/(3))
  }
  
  
  #ascending
  print("sim 1")
  #changes in beta
  count=0
  
  k1=function(t1,t2){k_exp(t1,t2,0.5,0.05)}
  c1=getCov(grid,k1)
  k2=function(t1,t2){k_exp(t1,t2,.75,0.05)}
  c2=getCov(grid,k2)
  k3=function(t1,t2){k_exp(t1,t2,1,0.05)}
  c3=getCov(grid,k3)
  k4=function(t1,t2){k_exp(t1,t2,1.25,0.05)}
  c4=getCov(grid,k4)
  k5=function(t1,t2){k_exp(t1,t2,1.5,0.05)}
  c5=getCov(grid,k5)
  k6=function(t1,t2){k_exp(t1,t2,1.75,0.05)}
  c6=getCov(grid,k6)
  
  
  # tmp=seq(0.5,1,length.out = 15)
  for(dist in c("N","t","sn")){
    fileName=paste0("FKWC_Results_cp/dist_",dist,"_mfn_depths_k_exp_Beta_asc_mult_N_",N,"_num_runs_200.Rda",sep="")
    runMVSim(N,c1,c2,c3,c4,c5,c6,grid,num_runs,fileName,dist=dist)
    count=count+1
    print(count/(3))
  }
  
  
  
  
  print("sim 2")
  #changes in alpha
  count=0
  k1=function(t1,t2){k_exp(t1,t2,0.5,.1)}
  c1=getCov(grid,k1)
  k2=function(t1,t2){k_exp(t1,t2,0.5,0.05)}
  c2=getCov(grid,k2)
  k3=function(t1,t2){k_exp(t1,t2,.5,.1)}
  c3=getCov(grid,k3)
  k4=function(t1,t2){k_exp(t1,t2,0.5,0.05)}
  c4=getCov(grid,k4)
  k5=function(t1,t2){k_exp(t1,t2,.5,.1)}
  c5=getCov(grid,k5)
  k6=function(t1,t2){k_exp(t1,t2,0.5,0.05)}
  c6=getCov(grid,k6)
  
  
  # tmp=seq(0.5,1,length.out = 15)
  for(dist in c("N","t","sn")){
    fileName=paste0("FKWC_Results_cp/dist_",dist,"_mfn_depths_k_exp_Alpha_alt_mult_N_",N,"_num_runs_200.Rda",sep="")
    runMVSim(N,c1,c2,c3,c4,c5,c6,grid,num_runs,fileName,dist=dist)
    count=count+1
    print(count/(3))
  }
  
  
  #ascending
  print("sim 2")
  #changes in beta
  count=0
  
  k1=function(t1,t2){k_exp(t1,t2,0.5,0.05)}
  c1=getCov(grid,k1)
  k2=function(t1,t2){k_exp(t1,t2,0.5,.075)}
  c2=getCov(grid,k2)
  k3=function(t1,t2){k_exp(t1,t2,0.5,.1)}
  c3=getCov(grid,k3)
  k4=function(t1,t2){k_exp(t1,t2,0.5,.125)}
  c4=getCov(grid,k4)
  k5=function(t1,t2){k_exp(t1,t2,0.5,.15)}
  c5=getCov(grid,k5)
  k6=function(t1,t2){k_exp(t1,t2,0.5,.175)}
  c6=getCov(grid,k6)
  
  
  # tmp=seq(0.5,1,length.out = 15)
  for(dist in c("N","t","sn")){
    fileName=paste0("FKWC_Results_cp/dist_",dist,"_mfn_depths_k_exp_Alpha_asc_mult_N_",N,"_num_runs_200.Rda",sep="")
    runMVSim(N,c1,c2,c3,c4,c5,c6,grid,num_runs,fileName,dist=dist)
    count=count+1
    print(count/(3))
  }
  
}









