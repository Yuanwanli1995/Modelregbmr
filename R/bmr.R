library(splines)
require(graphics)
library(quantreg)  
library(tcltk2)
library(lpme)
generatedata<- function(n,e,f){
  n<-n
  x<-sort(runif(n,0,1))
  if(e=="gauss"){
    error <- rnorm(n, mean= 0, sd = 1)
    mode<-rep(0,n)
  }else if(e=="mixgauss"){
    epsilon1<-rnorm(0.1*n,mean=-1,sd=2.5)
    epsilon2<-rnorm(0.9*n,mean=1,sd=0.5)
    epsilon<-c(epsilon1,epsilon2)
    error<-sample(epsilon,n,replace=TRUE)
    mode<-rep(1,n)
  }else if(e=="cauchy"){ 
    error <- rcauchy(n, location = 0, scale = .1)
    mode<-rep(0,n)
  }else if(e == "gamma"){
    error <- rgamma(n, shape=1,rate = 1)-0.3 
    mode<-rep(0,n)
  }else if(e == "beta"){
    error <- rbeta(n,5,3)
    mode<-rep(0,n)
  }
  if(f=="quodratic" ){
    func<-function(x){return(40*x^2-40*x+10)}
  }else if(f=="cubic"){
    func<-function(x){return(100*x^3-150*x^2+50*x+2)}
  }else if(f=="exp"){
    func<-function(x){return(x+exp(8*(x-0.5)^2))}
  }else if(f=="log"){
    func<-function(x){return(5*log(x,3)+10)}
  }else if(f=="sin"){
    func<-function(x){return(5*sin(20*x))}
  }else if(f=="mexicohat"){
    func<-function(x){return(-1+1.5*x+0.8*dnorm(x-0.6,0,0.04))}
  }else if(f=="steps"){
    func<-function(x) {
      ifelse(x < 0.25, 1, 
             ifelse(x < 0.5, 3,
                    ifelse(x < 0.75, 5,
                           ifelse(x <= 1, 1, NA))))
    }
  }
  y<- sapply(x, func)+error
  real<-sapply(x, func)+mode
  out <- list(x=x,y=y,real=real)
  return(out)
}

l2<-function(c){return(sqrt(t(c)%*%c))}#functions for EM 

phih<-function(i,X,y,beta_bmr,h){
  return(exp((y[i]-X[i,]%*%beta_bmr)^2/(-2*h^2))/(sqrt(2*pi)*h))}

ite_bmr<-function(x,X,y,beta_bmr,h){
  n<-length(x)
  W<-as.vector(phih(1:n,X=X,y=y,beta_bmr=beta_bmr,h=h))
  W<-diag(W)
  x.inv <- try(solve(t(X)%*%W%*%X),silent=T)
  if ('try-error' %in% class(x.inv)) return(beta_bmr)
  else return(solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%y)}

lmr<-function(x,y,beta_lmr,x0,z,h1,h2){
  pre<-beta_lmr[1]*z+beta_lmr[2]*(x-x0)+beta_lmr[3]*(x-x0)^2
  return(exp(-((x-x0)^2)/(2*h1^2)-((y-pre)^2)/(2*h2^2)))}

ite_lmr<-function(x,X,y,beta_lmr,x0,z,h1,h2){
  W<-diag(lmr(x,y,beta_lmr,x0,z,h1,h2))
  x.inv <- try(solve(t(X)%*%W%*%X),silent=T)
  if ('try-error' %in% class(x.inv)) return(beta_lmr)
  else return(solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%y)}

regression<-function(x,y,method,sp_1, sp_2){
  n<-length(x)
  if(method=="bmr"){
    h<-sp_1
    a<-sp_2
    aa=seq(0,n-1,a)
    knots1=c(0,0,0,0,x[aa],1,1,1,1) 
    X<-splineDesign(knots=knots1,x, outer.ok = F)
    beta_bmr<-rep(0,ncol(X))
    repeat{
      temp<-ite_bmr(x,X,y,beta_bmr,h)
      e=l2(temp-beta_bmr)
      beta_bmr<-temp
      if(e<0.001)break}
    fit_bmr<-X%*%beta_bmr
    fit<-fit_bmr
  }else if(method=="lmr"){
    #print("lmr calculating--")
    #print(strftime(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    fit_lmr<-vector()
    z<-rep(1,n)
    h1<-sp_1
    h2<-sp_2
    beta_lmr<-rep(1,3)
    for(k in 1:n){
      x0<-x[k]
      X<-cbind(z,x-x0,(x-x0)^2)
      repeat{
        temp<-ite_lmr(x,X,y,beta_lmr,x0,z,h1,h2)
        e=l2(temp-beta_lmr)
        # print(temp)
        # print(e)
        e2=l2(temp)
        beta_lmr<-temp
        if(e<0.01 | e2>100000)break}
      fit_lmr[k]<-beta_lmr[1]}
    fit<- fit_lmr
  }else if(method=="bqr"){
    fit<-fitted(rq(y~bs(x,df=3),tau=.5))
  }else if(method=="blse"){
    fit<-fitted(lm(y~bs(x,df=3)))
  }
  return(fit)
}

bandwidthselecte<-function(x,y,model,criterion,C=0.05,A=rev(seq(20,60,by=5)),H=rev(seq(0.01,2,by=0.01)),H1=-rev(seq(0.05,0.2,by=0.01)),H2=-rev(seq(0.01,1,by=0.02))){
  c<-(max(y)-min(y))*C
  if(model=="bmr"){
    S<<-matrix(rep(0,length(A)*length(H)),nrow = length(H))
    for(i in 1:length(H)){
      h<-H[i]
      for(j in 1:length(A)){
        a<-A[j]
        fit_bmr<-regression(x,y,method="bmr",h,a)
        if(criterion=="p" | criterion=="CV-mode" | criterion=="bootstrap"){
          S[i,j]<-sum(fit_bmr-c<y&y<fit_bmr+c)
        }else if(criterion=="mse"){
          S[i,j]<-sum((fit_bmr-y)^2)*(-1) 
        }
      }}
    max_index <- which.max(S)
    sp_1<-H[(max_index-1)%%length(H)+1]
    sp_2<-A[(max_index-1)%/%length(H)+1]
    if (criterion=="CV-mode"){
      hhxy=moderegbw(y, x, method="CV-mode", p.order=0,h1=H1, h2=H2)$bw;#get h2  #h1 x # h2 y 
      print(c("CV-mode数",hhxy))
      sp_1=hhxy[2]
    }else if(criterion=="bootstrap"){
      #hhxy=moderegbw(y, x, method="bootstrap", p.order=0,h1=H1, h2=H2,df=5, ncomp=5, nboot=5)$bw;#get h2
      hhxy=moderegbw(y, x, method="bootstrap",p.order=0,h1=H1, h2=H2)$bw;#get h2
      print(c("bootstrap",hhxy))
      sp_1=hhxy[2]
    }
  }else if(model=="lmr"){
    if (criterion=="CV-mode"){
      hhxy=moderegbw(y, x, method="CV-mode", p.order=0,h1=H1, h2=H2)$bw;
      sp_1=hhxy[1];sp_2=hhxy[2];
    }else if(criterion=="bootstrap"){
      hhxy=moderegbw(y, x, method="bootstrap", p.order=0, h1=H1, h2=H2)$bw;
      sp_1=hhxy[1];sp_2=hhxy[2];
    }else{
      S1<-matrix(rep(0,length(H1)*length(H2)),nrow = length(H1))
      for(i in 1:length(H1)){
        h1<-H1[i]
        for(j in 1:length(H2)){h2<-H2[j]
        fit_lmr<-regression(x,y,method="lmr",h1,h2)
        if(criterion=="p"){
          S1[i,j]<-sum(fit_lmr-c<y&y<fit_lmr+c)
        }else if(criterion=="mse"){
          S1[i,j]<-sum((fit_lmr-y)^2)*(-1)
        }
        }
      }
      max_index1 <- which.max(S1)
      sp_1<-H1[(max_index1-1)%%length(H1)+1]
      sp_2<-H2[(max_index1-1)%/%length(H1)+1]
    }
  }else if(model=="blse"){
    sp_1<--1;sp_2<--1
  }else if(model=="bqr"){
    sp_1<--2;sp_2<--2
  }
  out <- list(sp_1=sp_1,sp_2=sp_2)
  print(c("current paras",model,criterion,sp_1,sp_2,min(y),max(y)))
  return(out)
}

#Eeample
data<-generatedata(n=200,e="mixgauss", f="sin"); x<-data$x ; y<-data$y ; real <-data$real
best_spara<-bandwidthselecte(x,y,"bmr","p");sp_1<-best_spara$sp_1;sp_2<-best_spara$sp_2
fit<-regression(x,y,method="bmr",sp_1,sp_2)
plot(x,y, pch=1,cex.axis=1.5,cex.lab=1.5)
lines(x,fit,col="black",lwd=1,lty=3)
lines(x,real,col="red",lwd=0,lty=3)