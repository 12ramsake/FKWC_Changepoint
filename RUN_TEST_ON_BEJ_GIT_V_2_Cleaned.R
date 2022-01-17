#list of subjects
#
#sub08455
#sub08992
#sub08816
#sub34943
#sub12220
#sub06880
#

subjects=c('sub08455','sub08992','sub08816','sub34943','sub12220','sub06880')


#64x63x33 2 seconds apart
setwd("~/Functional Data Change Point Procedures/Beijing")

##PELT codes
Rcpp::sourceCpp('PELT_CPP.cpp')
#must load first
library(oro.nifti)
library(extrantsr)
library(neurobase)
library(ANTsR)
library(dplyr)
library(signal)
library(fda.usc)
library(fda)


subject=subjects[1]
for (subject in subjects){
  #extract image of proper subject
  #get list of images that could be fmri images
  sub_dir=paste("Beijing_",subject,sep="")
  pos_imgs=list.files(sub_dir,rec=T)
  sk_str=grepl("scan_rest",pos_imgs)
  image=pos_imgs[sk_str]
  t1 = neurobase::readnii(paste(sub_dir,"/",image,sep=""))
  
  
  
  #bias_correction, I think do this later
  # bc_t1 = bias_correct(file = t1, correction = "N4")
  
  
  curr_file=paste(sub_dir,"/",image,sep="")
  ##drop first 10 seconds
  
  ants_fmri = antsImageRead(curr_file)
  tr = 2 # 2 seconds
  first_scan = floor(10.0 / tr) + 1 # 10 seconds "stabilization of signal"
  sub_fmri = extrantsr::subset_4d(ants_fmri, first_scan:ntim(ants_fmri))
  

  
  ##prepare file for moco
  base_fname = paste(sub_dir,"/",sub_dir,sep="")
  avg_img = getAverageOfTimeSeries(sub_fmri)
  all_vox_dim = voxdim(sub_fmri )
  
  
  
  #####################
  # Motion Calculation
  ##################
  moco_file = paste0(base_fname, 
                     "_Motion_Params.rda")
  moco_fname = paste0(base_fname, "_moco_img.nii.gz")
  # sub_fmri=ants_fmri
  {
    if (all(file.exists(c(moco_file, 
                          moco_fname)))) { 
      load(moco_file)
      moco_img = antsImageRead(moco_fname)
      motion_res$moco_img = 
        moco_img
    } else {
      motion_res = 
        antsMotionCalculation(sub_fmri, 
                              fixed = avg_img, 
                              moreaccurate = 1,
                              txtype = "Rigid",
                              verbose = TRUE)
      save(motion_res, 
           file = moco_file)
      moco_img = 
        motion_res$moco_img
      antsImageWrite(moco_img, 
                     filename = moco_fname)
    }
  }
  
  #for plotting motion correction, assess later
  # moco_params =  motion_res$moco_params
  # # moco_params = moco_params %>% select(starts_with("MOCO"))
  # moco_params = moco_params[,3:8]
  # 
  # nuisanceVariables = moco_params
  # mp = round(moco_params, 4)
  # print(head(mp, 3))
  # 
  # 
  # 
  # mp = moco_params
  # mp[, 1:3] = mp[, 1:3] * 50
  # r = range(mp)
  # plot(mp[,1], type = "l", xlab = "Scan Number", main = "Motion Parameters",
  #      ylab = "Displacement (mm)",
  #      ylim = r * 1.25, 
  #      lwd = 2,
  #      cex.main = 2,
  #      cex.lab = 1.5,
  #      cex.axis = 1.25)
  # for (i in 2:ncol(mp)) {
  #   lines(mp[, i], col = i)
  # }
  # 
  
  
  #read in the motion corrected image
  #make an image mask, might not be useful
  moco_img = antsImageRead(moco_fname)
  moco_avg_img = getAverageOfTimeSeries(moco_img)
  maskImage = getMask(moco_avg_img, 
                      mean(moco_avg_img), 
                      Inf, cleanup = 2)
  mask_fname = paste0(base_fname, "_mask.nii.gz")
  antsImageWrite(maskImage, filename = mask_fname)
  
  mc=moco_img 
  
  #plot motion corrected image
  # neurobase::ortho2(mc)
  
  
  #run image through a high pass filter using signal package
  # 
  time=1:(dim(mc)[4])
  tr_mc=array(0,dim=dim(mc))
  bf <- butter(2, .1, type="high")
  high_pass_f=function(x){
    if(any(x>0)){
      val=signal::filtfilt(bf, x)
    }
    else
      val=rep(0,dim(tr_mc)[4])
    return(val)
  }
  # par(mfrow=c(1,2))
  # plot(mc[13,45,13,],type='l')
  # plot(high_pass_f(mc[13,45,13,]),type='l')
  # plot(high_pass_f(mc[sample(1:64,1),sample(1:64,1),sample(1:33,1),]),type='l')
  # 
  for( i in 1:64){
    for( j in 1:64){
      for( k in 1:33){
        tr_mc[i,j,k,]=high_pass_f(mc[i,j,k,])
      }
    }
  }
  # 
  # # turn back to nifiti image
  # tr_mc_i=neurobase::niftiarr(mc,tr_mc)
  # neurobase::ortho2(tr_mc_i)
  # 
  # mca=antsImageRead(mc)
  # tmp=ANTsR::getAverageOfTimeSeries(mca)
  # tmp2=tr_mc_i+tmp
  # plot(tr_mc_i[sample(1:64,1),sample(1:64,1),sample(1:33,1),],type='l')
  
  # 
  # mc2=mc
  # mc_filt=ANTsR::frequencyFilterfMRI(matrix(mc2[1,1,1,]),3,freqHi = 0.1)
  # 
  # tr_mc2=tr_mc
  # for( i in 1:64){
  #   for( j in 1:64){
  #     for( k in 1:33){
  #       tr_mc2[i,j,k,]=ANTsR::frequencyFilterfMRI(matrix(mc2[i,j,k,]),3,freqHi = 0.1)
  #     }
  #   }
  # }
  
  # tr_mc_2=neurobase::niftiarr(mc2,tr_mc2)
  # neurobase::ortho2(tr_mc_2)
  # neurobase::ortho2(mc)
  
  
  
  #we can try random, exp variogram directions
  #remove_edges
  trimm=16:(dim(tr_mc)[4]-15)
  tr_mc_no_edge=tr_mc[,,,trimm]

  # tr_mc_no_edge=tr_mc
  
  M=50

  
  #generating the proper unit functions
  
  dim_img=round(c(64,64,33))
  prod(dim_img)
  random_directions=array(0,dim=c(dim_img,M))
  
  #lets multiply together
  
  d1=fda.usc::rproc2fdata( M,t = seq(0,1,l=64),norm = TRUE)
  d2=fda.usc::rproc2fdata( M,t = seq(0,1,l=64),norm = TRUE)
  d3=fda.usc::rproc2fdata( M,t = seq(0,1,l=33),norm = TRUE)
  
  #Now combine into one function
  for(i in 1:64){
    for(j in 1:64){
      for(k in 1:33){
        random_directions[i,j,k,]=d1$data[,i]*d2$data[,j]*d3$data[,k]
      }
    }
  }
  
  #normalize
  for(i in 1:M){
    random_directions[,,,i]= random_directions[,,,i]/sqrt(sum(c(random_directions[,,,i])^2))
  }
  
  
  project=function(uv,img3d){
    num_bf=dim(uv)[1]
    proj=sum(uv*img3d)
    return(proj)
  }
  
  
  projections=sapply(1:M,function(rd){apply(tr_mc_no_edge,4,function(x){project(random_directions[,,,rd],x)})}); projections
  
  deriv_f_name=paste(base_fname,"_derivatives_HPfilter.rda")
  {
    if (file.exists(deriv_f_name)) {
      load(deriv_f_name)
    }
    else{
      # must compute derivative
      pic=function(x,pic_num=1){
        # print(x)
        xx=floor(x[1]*64)
        yy=floor(x[2]*64)
        zz=floor(x[3]*33)
        if(xx<1)
          xx=1
        if(zz<1)
          zz=1
        if(yy<1)
          yy=1
        
        if(xx>64)
          xx=64
        if(zz>33)
          zz=33
        if(yy>64)
          yy=64
        
        
        return(tr_mc_no_edge[xx,yy,zz,pic_num] )
        
      }
      
      all_pnts=expand.grid((1:64)/64,(1:64)/64,(1:33)/33)
      all_pnts2=expand.grid((1:64),(1:64),(1:33))
      
      get_derivatives=function(pic_num){
        
        deriv1=apply(all_pnts,1,function(z){numDeriv::grad(pic,unlist(z, use.names=FALSE),
                                                           method="Richardson",
                                                           method.args = list(eps=0.01,d=.1)
                                                           ,pic_num=pic_num)})
        
        dxx=dyy=dzz=tr_mc_no_edge[,,,pic_num]
        
        for(x in 1:nrow(all_pnts2)){
          y=unlist(all_pnts2[x,], use.names=FALSE);  
          dxx[y[1],y[2],y[3]] = deriv1[1,x]
          dyy[y[1],y[2],y[3]] = deriv1[2,x]
          dzz[y[1],y[2],y[3]] = deriv1[3,x]
        }
        
        return(list(dxx,dyy,dzz))
        
      }
      
      num_scans=dim(tr_mc_no_edge)[4]
      derivs=sapply(1:num_scans,get_derivatives)
      
      save(derivs,file=deriv_f_name)
      
    }
  }
  # 
  derivsx2=tr_mc_no_edge
  derivsy2=tr_mc_no_edge
  derivsz2=tr_mc_no_edge
  
  for(i in 1:(dim(tr_mc_no_edge)[4])){
    derivsx2[,,,i]=derivs[1,i][[1]]
    derivsy2[,,,i]=derivs[2,i][[1]]
    derivsz2[,,,i]=derivs[3,i][[1]]
  }
  
  projections_derivative_x=sapply(1:M,function(rd){apply(derivsx2,4,function(x){project(random_directions[,,,rd],x)})}); projections_derivative_x
  projections_derivative_y=sapply(1:M,function(rd){apply(derivsy2,4,function(x){project(random_directions[,,,rd],x)})}); projections_derivative_y
  projections_derivative_z=sapply(1:M,function(rd){apply(derivsz2,4,function(x){project(random_directions[,,,rd],x)})}); projections_derivative_z
  
  
  
  
  
  
  
  
  
  #for without derivatives
  # ranks_p=apply(projections,2,function(x){rank(x)})
  # ranks_p2=apply(projections,2,function(x){length(x)-rank(x)})
  # hds=mapply(min,ranks_p2,ranks_p)%>%matrix(nrow=dim(projections)[1],ncol=dim(projections)[2])
  # hdsm=apply(hds,1,mean)
  # # projections2=apply(projections,2,function(x){min(rank(x),225-rank(x))})
  # plot(hdsm/dim(projections)[1],type='l')
  # #Let's Get the Ranks
  # plot(rank(hdsm),type='l')
  
  
  
  #for with derivatives
  d_vals=sapply(1:M,function(x){ddalpha::depth.halfspace(cbind(projections[,x],projections_derivative_x[,x],projections_derivative_y[,x],projections_derivative_z[,x]),
                                                         cbind(projections[,x],projections_derivative_x[,x],projections_derivative_y[,x],projections_derivative_z[,x])) })
  
  
  
  # d_vals=sapply(1:M,function(x){ddalpha::depth.halfspace(cbind(projections[,x],projections_derivative[,x]),cbind(projections[,x],projections_derivative[,x])) })
  
  plot(rowMeans(d_vals),type='l')
  
  par(mfrow=c(1,1))
  plot(rank(rowMeans(d_vals)),type='l')
  # plot(rank(hdsm),type='l')
  
  
  
  
  cp=PELT_T(rank(rowMeans(d_vals)),.3*sqrt(201)+3.74)
  cp
  cp=cp[-1]
  
  # plot(rank(hdsm),type='l')
  # par(mfrow=c(1,1))
  plot(rank(rowMeans(d_vals)),type='l')
  abline(v=cp,col="red")
  
  # plot(tr_mc_no_edge[sample(1:64,1),sample(1:64,1),sample(1:33,1),],type='l')
  # abline(v=cp,col="red")
  rankz=rowMeans(d_vals)
  save(cp,rankz,file=paste(base_fname,"_changes.rda"))
}




subject=subjects[1] #
subject=subjects[2] #
subject=subjects[3] #
subject=subjects[4] #
subject=subjects[5] #
subject=subjects[6] #
for (subject in subjects){
  
  sub_dir=paste("Beijing_",subject,sep="")
  base_fname = paste(sub_dir,"/",sub_dir,sep="")
  load(paste(base_fname,"_changes.rda"))
  
  
  
  cp=PELT_T(rank(rankz),.4*sqrt(length(rank(rankz)))+3.74);  cp
  cp=cp[-1]
  fn=paste("CP_Graph",subject,".pdf",sep="")
  # Cairo::CairoPDF(file=fn,height=10,width=13)
  pdf(file=fn,height=10,width=13)
  plot(rank(rankz),type='l',lwd=3,main=subject,xlab="",bty="n",xaxt="n",yaxt="n",ylab="",cex.main=2)
  
  # abline(v=cp,col="magenta",lwd=3,lty=2)
  prev=1
  print(cp)
  cp2=c(cp,length(rankz)+1)
  for(i in 1:length(cp2)){
    lines(c(prev,(cp2[i]-1)),rep(mean(rank(rankz)[prev:(cp2[i]-1)]),2),col="magenta",lwd=4,lty=2)
    prev=cp2[i]
  }
  mtext("rank",2,line=1,cex=2)
  mtext("time",1,line=2,cex=2)
  axis(2,line=-1.5,cex.axis=2)
  axis(1,at=c(0,50,100,150,190),line=0,cex.axis=2)
  dev.off()
  # cp=PELT_T(rank(rankz)[15:220],.35*sqrt(length(rank(rankz)))+3.74);  cp
  # plot(rank(rankz),type='l')
  # abline(v=cp,col="red")
  
  # plot(rank(rankz),type='l')
  # abline(v=cp,col="red")
  
}

#AMOC model, run after


# compute test statistic and location of cp
run_test=function(rankz, boundary=0){
  n=length(rankz)
  Znk=function(kk){
    
    abs((((n)*(n^2-1)/12)^(-1/2))*sum(rankz[1:kk]-(n+1)/2))
    
  }
  
  
  
  Zns=sapply(boundary:(n-boundary),Znk)
  k1=which.max(Zns)
  k=(boundary:(n-boundary))[k1]
  Znt=Zns[k1]
  return(c(Znt,k))
}
#
#
#
# #discuss later
run_test(rank(rankz),0)

##BB cdf
ll=-1000:1000
BB_cdf<-function(q){sum(((-1)^ll)*exp(-2*ll^2*q^2))}


AMOC_results=NULL

for (subject in subjects){
  
  sub_dir=paste("Beijing_",subject,sep="")
  base_fname = paste(sub_dir,"/",sub_dir,sep="")
  load(paste(base_fname,"_changes.rda"))
  
  resu=run_test(rank(rankz))
  
  pval=1-BB_cdf(resu[1])
  
  print(resu)
  AMOC_results=rbind(AMOC_results,c(resu,pval))
  
  
  
}

AMOC_results
#CV is 1.35 at 95%



#epidemic model, run after
rankz=rank(rankz)

# compute test statistic and location of cp
run_test=function(rankz,boundary=10){
  n=length((rankz))
  sign=12/((n)*(n+1))
  mn=3*(n+1)
  Wk=function(k){
    k1=k[1]
    k2=k[2]
    # kruskal.test(rankz,g=c(rep(1,k1-1),rep(2,k2-k1),rep(1,n-k2+1)))$statistic
    t1=(sum(rankz[c(1:(k1-1),k2:n)])/sqrt(n-k2+k1))^2
    t2=(sum(rankz[k1:(k2-1)])/sqrt(k2-k1))^2
    return(sign*(t1+t2)-3*(n+1))
    # -sign*(cost_cpp(0,n-k2+k1-1, rankz[c(1:(k1-1),k2:n)])+cost_cpp(0,k2-k1-1, rankz[k1:(k2-1)]))-mn
    
    # return(sign*(  (n-k2+k1)*(mean(rankz[c(1:(k1-1),k2:n)])^2)+(k2-k1)*(mean(rankz[k1:(k2-1)])^2)   )-mn)
    
  }
  
  pairs=apply(combn(n,2),2,sort)
  pairs=pairs[,-(1:n)]
  pairs=pairs[,pairs[1,]!=n]
  pairs=pairs[,pairs[2,]!=n]
  
  #delete boundary points
  pairs=pairs[,pairs[1,]<(n-boundary)]
  pairs=pairs[,pairs[2,]<(n-boundary)]
  pairs=pairs[,pairs[1,]>(boundary+1)]
  pairs=pairs[,pairs[2,]>(boundary+1)]
  
  Zns=apply(pairs,2,Wk)
  ks=which.max(Zns)
  k=pairs[,ks]
  Znt=Zns[ks]
  return(c(Znt,k))
}
#
#
#
# #discuss later
#here we simulated the distribution of this test statistic, see the sim_bridge file 
load("~/Functional Data Change Point Procedures/Simulation Study/Epidemic/epid_dist2.rda")

run_test(rank(rankz))
EP_results=NULL
for (subject in subjects){
  
  sub_dir=paste("Beijing_",subject,sep="")
  base_fname = paste(sub_dir,"/",sub_dir,sep="")
  load(paste(base_fname,"_changes.rda"))
  
  resu=run_test(rank(rankz))
  
  pval=mean(maxes>resu[1])
  
  print(resu)
  EP_results=rbind(EP_results,c(resu,pval))
}


EP_results


quantile(maxes,.95)


#CV is 20.5 at 95%



#Make a table
AMOC_results[,3]=round(AMOC_results[,3],4)
tog=cbind(AMOC_results,EP_results)
tog=tog[,-c(1,4)]
colnames(tog)=c('Estimate','p-value','Estimate 1','Estimate 2','p-value')
rownames(tog)=subjects
xtable::xtable(tog)








# [1] 11.31456 34.00000 54.00000    # has outlier at beginning, hypothesis is not rejected
# [1]  61.19692  36.00000 186.00000 # we see a change;  but they don't
# [1]  11.59917  37.00000 186.00000 # we see no change ; stoehr found epidemic
# [1]  23.45679  21.00000 173.00000 # we do not see a change, they do not
# [1]  16.56529  31.00000 127.00000 # we see no change; they see change
# [1]  42.645  24.000 117.000       # we see a change; they see change
















