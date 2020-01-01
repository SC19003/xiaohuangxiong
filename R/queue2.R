#' @title Use four inputs to predict response using R.
#' @description The prediction model is described in http://www.babelgraph.org/wp/?p=358.
#' @param lambda the successive arrival intervals of customers to the system are independent and subject to an exponential distribution with a parameter of $\lambda$. (numeric)
#' @param mu The service time of the service desk is also iid., and obeys the exponential distribution of parameter $\mu$. (numeric)
#' @param T System service time (numeric)
#' @param S the number of the service desks (numeric)
#' @return  the length Ls, the average stay time Ws and the probability of customer waiting Pwait.
#' @examples
#' \dontrun{
#' lambda <- 4
#' mu <- 6
#' S <- 4
#' res <- queue2(4,6,10000,4)
#' }
#' @export
queue2<-function(lambda, mu, T, S){
  k<-0; wt<-0; wn<-0; ws<-0
  tp<-0; nA<-0; t<-0
  r<-runif(1); tA<--1/lambda*log(r)
  tD<-rep(Inf, S); SS<-rep(0, S+1)
  repeat{
    t1<-if(SS[1]==0) Inf else min(tD)
    i1<-if(SS[1]==0) 1 else which.min(tD)
    k<-k+1; wt[k]<-t; wn[k]<-SS[1]
    if (tA < T){
      ws[k]<-min(tA, t1)-t
      if (tA < t1){
        t<-tA; nA<-nA+1
        r<-runif(1); tA<-t-1/lambda*log(r)
        n<-SS[1]; SS[1]<-n+1
        for (i in 1:S){
          if (SS[1+i]==0){
            SS[1+i]<-1
            r<-runif(1); tD[i]<-t-1/mu*log(r)
            break
          }
        }
        
      }else{
        t<-t1; n<-SS[1]; SS[1]<-n-1
        if (n==1){
          SS[2:(S+1)]<-0; tD[1:S]<-Inf
        }else if (n<=S){
          SS[1+i1]<-0; tD[i1]<-Inf
        }else{
          r<-runif(1); tD[i1]<-t-1/mu*log(r)
        } }
    }else{
      ws[k]<- if( t1==Inf) 0 else t1-t
      n<-SS[1]
      if (n>0){
        t<-t1; SS[1]<-n-1;
        if (n==1){
          SS[2:(S+1)]<-0; tD[1:S]<-Inf
        }else if (n<=S){
          SS[1+i1]<-0; tD[i1]<-Inf
        }else{
          r<-runif(1); tD[i1]<-t-1/mu*log(r)
        }
      }else
        tp<-1
    }
    if (tp==1) break
  }
  data.frame(Ls=sum(ws*wn)/t, Ws=sum(ws*wn)/nA,
             Pwait=sum(ws[wn>=S])/t)
}
