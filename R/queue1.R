#' @title Use three inputs to predict response using R.
#' @description The prediction model is described in http://www.babelgraph.org/wp/?p=358.
#' @param lambda the successive arrival intervals of customers to the system are independent and subject to an exponential distribution with a parameter of $\lambda$. (numeric)
#' @param mu The service time of the service desk is also iid., and obeys the exponential distribution of parameter $\mu$. (numeric)
#' @param T System service time (numeric)
#' @return  the length Ls, the average stay time Ws and the probability of customer waiting Pwait.
#' @examples
#' \dontrun{
#' lambda <- 4
#' mu <- 6
#' T <- 1000
#' res <- queue1(4,6,1000)
#' }
#' @export
queue1<-function(lambda, mu, T){
  k<-0; wt<-0; wn<-0; ws<-0;
  tp<-0; nA<-0; n<-0; t<-0
  r<-runif(1); tA<--1/lambda*log(r); tD<-Inf
  repeat{
    k<-k+1; wt[k]<-t; wn[k]<-n
    if (tA < T){
      ws[k]<-min(tA, tD)-t
      if (tA < tD){
        t<-tA; n<-n+1; nA<-nA+1
        r<-runif(1); tA<-t-1/lambda*log(r)
        if (n==1){
          r<-runif(1); tD<-t-1/mu*log(r)
        }
      }else{
        t<-tD; n<-n-1
        if (n==0){
          tD<-Inf
        }else{
          r<-runif(1); tD<-t-1/mu*log(r)
        } }
    }else{
      ws[k]<-if(tD==Inf) 0 else tD-t
      if (n>0){
        t<-tD; n<-n-1
        if (n>0){
          r<-runif(1); tD<-t-1/mu*log(r)
        }
      }else
        tp<-1
    }
    
    if (tp==1) break
  }
  data.frame(Ls=sum(ws*wn)/t, Ws=sum(ws*wn)/nA,
             Pwait=sum(ws[wn>=1])/t)
}
