# M-quantile Spectral Discriminant Analysis (`mqper` and related functions)

A collection of R functions for performing robust discriminant analysis in the frequency domain using M-quantile regression methods.

---

## Description
This suite of functions implements M-quantile spectral analysis and discriminant classification for time series data. The main components are:

1. **`mqper`**: Computes M-quantile periodograms for robust spectral estimation
2. **`lperd.mtm`**: Smooths and logs the periodogram
3. **`cep.mtm`**: Computes cepstral coefficients from the log-periodogram
4. **`cep.lda`**: Main function performing discriminant analysis using cepstral features

Key features:
- Robust spectral estimation using M-quantile regression
- Cepstral coefficient extraction for feature representation
- Linear discriminant analysis (LDA) for classification
- Automatic selection of optimal number of cepstral coefficients
- Cross-validation support

---

## Installation
Ensure the following package is installed:
```R
install.packages("MASS")

# Generate sample time series data
set.seed(123)
# Download the following functions from \url{https://github.com/fppatrick}
source("mquantile.R") # to estimate the M-quantile regression model
source("mqper.R") # to estimate the M-quantile periodogram
source("mquantile.rlda.R") 
# Function to cfreate a data frame with several univariate time series
r.cond.ar2 <- function(N, nj=1, r.phi1, r.phi2, r.sig2){
        X = ts(matrix(0,ncol=nj,nrow=N))
        parms = matrix(0,ncol=nj,nrow=3)
        for(j in 1:nj){
                phi1 = runif(1,min=min(r.phi1), max=max(r.phi1) )
                phi2 = runif(1,min=min(r.phi2), max=max(r.phi2) )
                sigma2 = runif(1,min=min(r.sig2), max=max(r.sig2) )
                Xj = arima.sim(list(order=c(2,0,0),ar=c(phi1,phi2)),sd=sqrt(sigma2),n=N)
                X[,j] = Xj
                parms[,j] = c(phi1,phi2,sigma2)
        }
        z=list(X=X, parms=parms)
}
# number of series in training data
nj = 50
# length of time series 
N = 250 
# Define the train data sets
traindata1 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.01,.7),r.phi2=c(-.12,-.06),r.sig2=c(.1,10))
traindata2 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.5,1.2),r.phi2=c(-.36,-.25),r.sig2=c(.1,10))
traindata3 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.9,1.5),r.phi2=c(-.56,-.75),r.sig2=c(.1,10))
train <- cbind(traindata1$X,traindata2$X,traindata3$X)
# Create the labels
y <- c(rep(1,nj),rep(2,nj),rep(3,nj))
# Fit the discriminant analysis
fitrob <- cep.lda(y,train,tau = 0.5)
# Calculate the classifictaion rate
mean(fit$predict$class == y)
# Create the confusion matrix
table(y,fit$predict$class)
