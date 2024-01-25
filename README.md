# Modelregbmr
A nonparametric model based on B-splines is given for modal regression.


## Documentation
B-splines Mode Regression and Cross-Validation prediction  Bandwidth Selector for B-splines Mode Regression.

## Usage
#### Eeample
data<-generatedata(400,"mixgauss", "log"); x<-data$x ; y<-data$y ; real <-data$real

best_spara<-bandwidthselecte(x,y,"bmr","p");sp_1<-best_spara$sp_1;sp_2<-best_spara$sp_2

fit<-regression(x,y,method="bmr",sp_1,sp_2)

fit<-regression(x,y,method="lmr",sp_1,sp_2)

plot(x,y, pch=1,cex.axis=1.5,cex.lab=1.5)

lines(x,fit,col="black",lwd=1,lty=3)

lines(x,real,col="red",lwd=0,lty=3)


##  Arguments
**x**
an n by 1 response vector.

**y**
an n by 1 predictor vector.

**method**	
method="bmr" is for B-splines hyperparameter k in Mode Regression; method="lmr" is for Local polynomials Mode Regression;  For non-measurement error models, method="bmr" is recommended.

**C**
Percentage of abs((max(y)-min(y))),the constant hyperparameter of CVp

**A**	
bandwidth vector for B-splines hyperparameter k in Mode Regression ; default is NULL, and k is chosen automatically.  It is recommended to carefully specify a fine grid for k.

**H**
bandwidth vector for h in  B-splines modal regression; default is NULL, and h is chosen automatically.  It is recommended to carefully specify a fine grid for h.

**H1**	
bandwidth vector for h1 in LPM; default is NULL, and h1 is chosen automatically.  It is recommended to carefully specify a fine grid for h1.

**H2**
bandwidth vector for h2  in LPM; default is NULL, and h2 is chosen automatically.  It is recommended to carefully specify a fine grid for h2.

**n**	
length of vector x or y

**e** ="gauss","mixgauss","cauchy", "gamma","beta"  is for diffrent noise; 

**f** ="quodratic","cubic","exp", "log","sin","mexicohat","steps" is for diffrent simulation function; 

## Author(s)

Lianqiang Yang, Wanli Yuan, Shijie Wang


