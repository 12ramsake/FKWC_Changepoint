# This script is for the main Epidemic simulation in the paper.
#It only computes the depth values for each of the 200 runs . 
#You must then run the Epidemic or PELT algorithm on the ranks of the depth values after




# 3at least 10% apart
compute_location_change_points=function(N){
  min_length=N*.1
  candidate=sort(runif(2,.1,.9))
  while(diff(candidate)<.1)
    candidate=sort(runif(2,.1,.9))
  changes=floor(candidate*N)
  N1=changes[1]
  N2=changes[2]-changes[1]
  N3=N-changes[2]
  return(c(N1,N2,N3))
}


set.seed(440)
changes_100=replicate(num_runs, compute_location_change_points(100))
changes_200=replicate(num_runs, compute_location_change_points(200))
changes_500=replicate(num_runs, compute_location_change_points(500))

#now remove the seed!
set.seed(NULL)


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
#generate data from spec. covariance kernel, under epidemic model
gfd_k_3<-function(sam_size,
                  res=1e2, 
                  Cov1,
                  Cov2,
                  dist="t",
                  del=0.9){
  
  
  grid = seq( 0, 1, length.out = res )
  
  #compute the placement of the epidemic model
  #
  # locations=compute_location_change_points(N)
  N1=sam_size[1]
  N2=sam_size[2]
  N3=sam_size[3]
  
  if(dist=="N"){
    Data1 = generate_gauss_fdata( N1,  rep(0,res),Cov = Cov1 )
    Data2 = generate_gauss_fdata( N2, rep(0,res), Cov = Cov2 )
    Data3 = generate_gauss_fdata( N3, rep(0,res), Cov = Cov1 )
  }
  else if(dist=="t"){
    Data1 = rmvt( N1,  sigma =  Cov1/3,df=3 )
    Data2 = rmvt( N2,  sigma =  Cov2/3,df=3 )
    Data3 = rmvt( N3,  sigma =  Cov1/3 ,df=3)
  }
  else{
    # Data1 = rdmsn(N1, res, rep(0,res),Cov1, del=rep(del,res))
    # Data2 = rdmsn(N2, res, rep(0,res),Cov2, del=rep(del,res))
    # Data3 = rdmsn(N3, res, rep(0,res),Cov1, del=rep(del,res))
    
    
    cp1=list(mean=rep(0,nrow(Cov1)), var.cov=Cov1, gamma1=rep(del,nrow(Cov1))/nrow(Cov1))
    cp2=list(mean=rep(0,nrow(Cov2)), var.cov=Cov2, gamma1=rep(del,nrow(Cov2))/nrow(Cov2))
    cp3=list(mean=rep(0,nrow(Cov1)), var.cov=Cov1, gamma1=rep(del,nrow(Cov1))/nrow(Cov1))
    dp1=cp2dp(cp1, "SN")
    dp2=cp2dp(cp2, "SN")
    dp3=cp2dp(cp3, "SN")
    
    Data1 = rmsn(N1, dp=dp1)
    Data2 = rmsn(N2, dp=dp2)
    Data3 = rmsn(N3, dp=dp3)
    
  }
  
  return(list(argvals =  grid, mdata=rbind(Data1,Data2,Data3)))
  
}



#simulate epidemic model and then compute depth functions
simulate_mv_depths<-function(sam_size,c1,c2,res=100,
                             grid=seq( 0, 1, length.out = res ),
                             dist="t",
                             del=0.9
){
  
  dat=gfd_k_3(sam_size,Cov1=c1,Cov2=c2,res=res,dist=dist,del=del)
  rownames(dat$mdata)=as.character(1:sum(sam_size))
  fdat<-fdata(dat$mdata,dat$argvals)
  
  derivs<-derivcurves(dat$mdata)
  fderivs<-fdata(derivs,dat$argvals)
  
  # depth.RPD() implements a depth measure based on random projections possibly using several derivatives (see Cuevas et al. 2007).
  RPD_depth=RPD(dat$mdata)
  #fraiman and muniz, really slow
  FM_depth=depth.FM(fdat)$dep
  LTR_depth=c(norm.fdata(fdat))
  
  
  RPD_depth_d=RPDd(dat$mdata,derivs)
  #fraiman and muniz, really slow
  FM_depth_d=FMp(dat$mdata,derivs)
  
  #LTR
  LTR_depth_d=c(norm.fdata(fdat))+c(norm.fdata(fderivs))
  
  
  depths=data.frame(FM_depth,RPD_depth,LTR_depth,FM_depth_d,RPD_depth_d,LTR_depth_d)
  
  return(depths)
  
}



# runMVSim(N,c1,c2,grid,num_runs,fileName,dist=dist)
runMVSim<-function(N=50,
                   Cov1,
                   Cov2,
                   grid,num_runs,
                   fileName,
                   dist="t",
                   del=0.9){
  
  
  
  no_cores<-detectCores()-1
  
  cl <- makeCluster(no_cores,type="FORK")   
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  
  if(N==100)
    depth_values=try({foreach(i=1:num_runs) %dopar% {simulate_mv_depths(changes_100[,i],Cov1,Cov2,res=1e2,del=del,dist=dist,grid=grid)}})
  if(N==200)
    depth_values=try({foreach(i=1:num_runs) %dopar% {simulate_mv_depths(changes_200[,i],Cov1,Cov2,res=1e2,del=del,dist=dist,grid=grid)}})
  if(N==500)
    depth_values=try({foreach(i=1:num_runs) %dopar% {simulate_mv_depths(changes_500[,i],Cov1,Cov2,res=1e2,del=del,dist=dist,grid=grid)}})
  
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





# N

for(N in c(100,200,500)){
  # N2=N1
  # N3=N1
  print("sim 1")
  #changes in beta
  count=0
  #outside epidemic
  k1=function(t1,t2){k_exp(t1,t2,0.5,0.05)}
  c1=getCov(grid,k1)
  
  #inside epidemic
  k2=function(t1,t2){k_exp(t1,t2,1,0.05)}
  c2=getCov(grid,k2)
  #
  for(dist in c('N','t','sn')){
  
    fileName=paste0("FKWC_Results_cp/dist_",dist,"_mfn_depths_k_exp_Beta_Epidemic_N_",N,"_num_runs_200.Rda",sep="")
    runMVSim(N,c1,c2,grid,num_runs,fileName,dist="t")
    count=count+1
    print(count/(3))
  }
  
  
  
  
  print("sim 2")
  #changes in alpha
  count=0
  k1=function(t1,t2){k_exp(t1,t2,0.5,0.05)}
  c1=getCov(grid,k1)
  
  k2=function(t1,t2){k_exp(t1,t2,0.5,0.1)}
  c2=getCov(grid,k2)
  
  for(dist in c('N','t','sn')){
    fileName=paste0("FKWC_Results_cp/dist_",dist,"_mfn_depths_k_exp_Alpha_Epidemic_N_",N,"_num_runs_200.Rda",sep="")      
    runMVSim(N,c1,c2,grid,num_runs,fileName,dist="t")
    count=count+1
    print(count/(3))
  }
  
  
  
  #####under null
  
  
  
  print("sim 3")
  #changes in alpha
  count=0
  k1=function(t1,t2){k_exp(t1,t2,0.5,0.05)}
  c1=getCov(grid,k1)
  
  k2=function(t1,t2){k_exp(t1,t2,0.5,0.05)}
  c2=getCov(grid,k2)
  
  for(dist in c('N','t','sn')){
    fileName=paste0("FKWC_Results_cp/dist_",dist,"_mfn_depths_k_exp_NULL_Epidemic_N_",N,"_num_runs_200.Rda",sep="")      
    runMVSim(N,c1,c2,grid,num_runs,fileName,dist="t")
    count=count+1
    print(count/(3))
    
  }
  
  
  
}