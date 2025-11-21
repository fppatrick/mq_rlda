#' M-Quantile spectral discriminant analysis
#' Description: This function fits a robust discriminant analysis in the
#'  domain using the M-quantile regression method.
#'
#' @author: Patrick Ferreira Patrocinio
#' parameters:
#' y: Class labels vector
#' @param x: Time series matrix
#' @param tau: Quantile level (0 < tau < 1)
#' @param xNew: Optional new data for prediction
#' @param L: Number of cepstral coefficients (if FALSE, uses CV)
#' @param mcep: Maximum number of cepstral coefficients to consider
#' @param cv: Whether to perform cross-validation
#' @param tol: Tolerance for LDA

mqper <- function(series, tau) {
  n <- length(series)
  g <- n %/% 2
  per <- FFT <- rep(0,g)
  for (j in 1:(g-1)) {
    X1 <- X2 <- NULL
    w <- 2 * pi * j / n
    X1 <- cos(w * 1:n)
    X2 <- sin(w * 1:n)
    MX <- cbind(X1, X2)
    fitrob <- mqlm(MX, series, maxit=100, q=tau, k=1.345)
    FFT[j] <- sqrt(n / (8 * pi)) * complex(real = fitrob$coefficients[2], imaginary = -fitrob$coefficients[3])
    
  }
  w <- 2 * pi * seq_len(g) / n
  X1 <- cos(w * 1:n)
  fitrob <- mqlm(MX, series, maxit=100, q=tau, k=1.345)
  FFT[g] <- sqrt(n / (8 * pi)) * complex(real = fitrob$coefficients[2], imaginary = -0)
  per <- Mod(FFT)^2
  
  perior <- c(per[-1], rev(per))
  
  w <- 2 * pi * seq.int(g) / n
  return(list(spec = perior, freq = w))
}


lperd.mtm <-function(x, tau){
        mtm = mqper(x, tau)
        spec = smooth.spline(mtm$spec)$yin
        lspec = log(spec)
        fr = mtm$freq
        z=list(freq=fr, lspec=lspec, spec=spec)
}


cep.mtm <-function(x, tau){
  n = length(x)
  lpa = lperd.mtm(x, tau)
  lp1 = lpa$lspec
  cp = fft(c(lp1[-n], 0))/(2 * pi)
  z=list(quef=0:(n-1), cep=Re(cp), freq=lpa$freq, lspec=lpa$lspec)
  z
}

cep.get <- function(y, x, tau){
  if(dim(x)[2] != length(y))
    stop("\n Number of time series and group information (y) must be the same \n")
  n <- dim(x)[2]
  N <- dim(x)[1]
  log.spec.hat <- matrix(0,nrow=n,ncol=N)
  cep.hat <- matrix(0,nrow=n,ncol=N)

  results <- apply(x, 2, function(col) cep.mtm(col, tau))
  log.spec.hat <- do.call(rbind, lapply(results, `[[`, "lspec"))
  cep.hat <- do.call(rbind, lapply(results, `[[`, "cep"))
  D.hat <- cbind(data.frame(cep.hat),y)
  colnames(D.hat)[1:N]= paste("C", 0:(N-1), sep="")
  return(D.hat)
}


# Cross-validation
Lopt.get <- function(data, mcep=mcep){
  cvK <- array(0,dim=mcep)
  for(k in 1:mcep){
    b <- as.formula(paste("y ~ ",paste(colnames(data[,1:(k+1)]), collapse="+"),sep = ""))
    C.lda.pred <- lda(b , data=data, CV=TRUE)
    cvK[k] <- mean(C.lda.pred$class==data$y)
  }
  Lopt <- min(which(cvK == max(cvK)))
  return(Lopt)
}


cep.lda <- function(y,x,tau,xNew=NULL,L=FALSE,mcep=10,cv=FALSE,tol=1.0e-4
){
  if(is.matrix(x)==FALSE)
    stop("\n x must be a matrix or cannot be a signle time series")
  if(ncol(x) < 2)
    stop("\n Must have at least two time series observations")
  if(length(y) != dim(x)[2])
    stop("\n Length of y and numbers of time series in x must agree")
  if(is.null(xNew) == 0 & cv != 0)
    stop("\n Cannot do precition of new data and leave-out-one prediction of training data.")
  
  D.hat0 <- cep.get(y,x, tau)
  
  if(isTRUE(L)){
    Lopt <- L
  }else{
    Lopt <- Lopt.get(D.hat0,mcep)
  }
  
  b <- as.formula(paste("y ~ ",paste(colnames(D.hat0[,1:(Lopt+1)]), collapse="+"),sep = ""))
  if(isTRUE(cv)){
    C.lda <- lda(b, data=D.hat0, CV=TRUE,tol=tol)
    pre = list(class=C.lda$class, posterior=C.lda$posterior)
    cep.lda <- list(C.lda=C.lda, cep.data=D.hat0, Lopt=Lopt, lspec = NULL, predict=pre)
  } else{
    C.lda <- lda(b, data=D.hat0, CV=FALSE,tol=tol)

    Q = min(Lopt,length(unique(y))-1)
    frq <- seq(from=0, to=.5, by=1/(dim(D.hat0)[2]-1))
    dsc <- matrix(NA, nrow=length(frq), ncol=min(Lopt,length(unique(y))-1))
    for(q in 1:Q){
      dsc[,q] = recon(C.lda$scaling[,q],frq)
    }
    lspec <- list(dsc=dsc, frq=frq)               

    if(is.null(xNew)==0){
      pre <- predict(C.lda, cep.get(array(0,dim=dim(xNew)[2]),xNew), prior=C.lda$prior)
    }else{
      pre <- predict(C.lda, D.hat0, prior=C.lda$prior)
    }
    cep.lda <- list(C.lda=C.lda, cep.data=D.hat0, Lopt=Lopt, lspec=lspec, predict=pre)
  }
  class(cep.lda) <- "ceplda"
  return(cep.lda)
}



recon <- function(wts, fqs){
  nw = length(wts)
  FBmat = matrix(0,nrow=length(fqs),ncol=nw[1])
  FBmat[,1] = rep(1,length(fqs))
  for(j in 2:nw[1]){
    FBmat[,j] = cos(2*pi*(j-1)*fqs)
  }
  dsc = FBmat %*% wts
  dsc
}
