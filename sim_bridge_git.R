# install.packages('e1071')

gridd=1000
bridges=replicate(10000,e1071::rbridge(1,gridd))
dim(bridges)
beepr::beep()
install.packages('sde')

bridges=replicate(10000,sde::BBridge(x=0, y=0, t0=0, T=1, N=gridd))
dim(bridges)

maximize_bb=function(bi){
  print(bi)
  bridge_fn=Vectorize(function(k1,k2,bi=1){
    
    ((1-k2/gridd+k1/gridd)^(-1)+(k2/gridd-k1/gridd)^(-1))*(bridges[k2,bi]-bridges[k1,bi])^2
    
  },vectorize.args = c('k1','k2'))
  
  tst=outer(2:gridd,2:gridd,bridge_fn,bi=bi)
  tst[lower.tri(tst,T)]=-gridd
  tst2=unlist(tst)
  tst2[is.nan(tst2)]=-gridd
  # max(tst2)
  val=max(tst2)
  print(val)
  return(val)
}

maxes=sapply(1:10000,maximize_bb)
quantile(maxes,.95)
# 20.49254 
# 20.4


dirr<-"~/Functional Data Change Point Procedures/Simulation Study/Epidemic/"
setwd(dirr)

quantile(maxes,.95)
save(maxes,file="epid_dist2.rda")







