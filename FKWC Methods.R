#given the depth ranks, the FKWC methods are below

#you may need these

library(MASS)
library(stringi)
library(xtable)
library(rlist)
library(latex2exp)



####################### AMOC #######################

#estimate cdf of BB cdf

ll=-1000:1000
cdf<-function(q){sum(((-1)^ll)*exp(-2*ll^2*q^2))}
qBB<-function(p){return(uniroot(function(q){cdf(q)-p},c(1,1.5))$root)}
qBB(0.95)
thresh=qBB(0.95)
thresh


#runs the AMOC test on a set of ranks
boundary=1
run_test=function(rankz){
  
  Znk=function(kk){
    
    abs((((n)*(n^2-1)/12)^(-1/2))*sum(rankz[1:kk]-(n+1)/2))
    
  }
  
  Zns=sapply(boundary:(n-boundary),Znk)
  k1=which.max(Zns)
  k=(boundary:(n-boundary))[k1]
  Znt=Zns[k1]
  return(c(Znt,k))
}
# Znt>thresh?

####################### Epidemic #######################

#estimate cdf

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
# 20.4



#my C++ implementation of KW PELT should be sourced
Rcpp::sourceCpp('PELT_CPP.cpp')

#compute test statistic and location of cp
run_test=function(rankz){
  n=length(rankz)
  sign=12/((n)*(n+1))
  mn=3*(n+1)
  
  Wk=function(k){
    k1=k[1]
    k2=k[2]
    
    -sign*(cost_cpp(0,n-k2+k1-1, rankz[c(1:(k1-1),k2:n)])+cost_cpp(0,k2-k1-1, rankz[k1:(k2-1)]))-mn
    

  }
  
  pairs=apply(combn(n,2),2,sort)
  
  Zns=apply(pairs,2,Wk)
  ks=which.max(Zns)
  k=pairs[,ks]
  Znt=Zns[ks]
  return(c(Znt,k))
}


# Znt>thresh?

####################### Multiple change-point #######################

#my C++ implementation of KW PELT should be sourced
Rcpp::sourceCpp('PELT_CPP.cpp')

#Then use PELT_T to run the algorithm
h=0.3
lambdan=3.74+sqrt(N)*h
PELT_T(ranks,beta=lambdan)

