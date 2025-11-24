#'@title M-quantile spectral discriminant analysis
#'
#' Helper functions to compute M-quantile periodograms, cepstral
#' coefficients and perform cepstral LDA (spectral discriminant analysis).
#'
#' @author: Patrick Ferreira Patrocinio
#' Note: These functions depend on a working `mqper()` and `mqlm()` implemented
#'       elsewhere in your package / workspace.
#'
#' Compute (smoothed) log-periodogram from mqper output
#'
#' This wrapper calls \code{mqper(x, tau)} and smooths the raw periodogram
#' across frequency using \code{stats::smooth.spline}. It returns the
#' smoothed spectrum, its log and the original frequency vector.
#'
#' @param x Numeric vector: univariate time series.
#' @param tau Numeric in (0,1): M-quantile level.
#' @return A list with elements \code{freq}, \code{spec} (smoothed spectrum) and \code{lspec} (log of smoothed spectrum).
#' @examples
#' \dontrun{ lperd.mtm(rnorm(128), 0.5) }
lperd.mtm <- function(x, tau) {
  if (!is.numeric(x)) {
    stop("x must be a numeric vector")
  }
  if (!is.numeric(tau) || tau <= 0 || tau >= 1) {
    stop("tau must be in (0,1)")
  }

  # mqper is expected to return a list with at least 'perior' (periodogram) and 'freq'
  mtm <- mqper(x, tau)

  # --- Safety: Accept both 'perior' and 'spec' naming conventions if present
  if (!is.null(mtm$perior)) {
    raw_spec <- as.numeric(mtm$perior)
  } else {
    stop("mqper() output must contain 'perior' periodogram values")
  }

  freq <- as.numeric(mtm$freq)
  if (length(freq) != length(raw_spec)) {
    stop("length of frequencies and periodogram do not match")
  }

  # Smooth the periodogram across frequency (predict at original freq locations)
  ss <- stats::smooth.spline(x = freq, y = raw_spec)
  #spec_smoothed <- stats::predict(ss, freq)$y

  # Avoid non-positive values before log
  ss$yin[ss$yin <= 0] <- .Machine$double.eps

  lspec <- log(ss$yin)

  list(freq = freq, spec = ss$yin, lspec = lspec)
}

#' Cepstral coefficients from (smoothed) log-periodogram
#'
#' Computes cepstral coefficients as the real part of the FFT of the
#' (shifted) log-spectrum. For a time series of length N we compute the
#' FFT of [lspec_1, ..., lspec_{N-1}, 0] and return the real part.
#'
#' @param x Numeric vector: univariate time series of length N.
#' @param tau Numeric in (0,1): M-quantile level.
#' @return A list containing \code{quef} (0:(N-1)), \code{cep} (real cepstral coefficients),
#'         \code{freq} and \code{lspec}.
#' @examples
#' \dontrun{ cep.mtm(rnorm(128), 0.5) }
cep.mtm <- function(x, tau) {
  if (!is.numeric(x)) stop("x must be numeric")
  N <- length(x)

  lpa <- lperd.mtm(x, tau)
  lp1 <- as.numeric(lpa$lspec)

  # Form vector for FFT: drop last element and append 0 (as original code did)
  fft_input <- c(lp1[-N], 0)
  cp <- stats::fft(fft_input) / (2 * pi)

  list(quef = 0:(N - 1L),
       cep  = Re(cp),
       freq = lpa$freq,
       lspec = lpa$lspec)
}

#' Build cepstral data frame for multiple time series
#'
#' Given a matrix \code{x} whose columns are individual time series (each of length N),
#' and a class label vector \code{y} of the same length as the number of columns,
#' returns a data frame with cepstral coefficients (C0..C_{N-1}) and the class label.
#'
#' @param y Factor or vector of class labels (length equals number of columns of x).
#' @param x Numeric matrix with time along rows and series across columns (N x n_series).
#' @param tau Numeric: M-quantile level.
#' @return A data.frame with columns C0..C_{N-1} and column y.
#' @examples
#' \dontrun{
#'   X <- matrix(rnorm(128*10), nrow = 128, ncol = 10)
#'   cep.get(rep(1:2, each=5), X, 0.5)
#' }
cep.get <- function(y, x, tau) {
  if (!is.matrix(x)) stop("x must be a matrix with time along rows and series across columns")
  n_series <- ncol(x)
  if (length(y) != n_series) stop("Number of class labels (y) must equal number of columns in x")

  # Apply cepstral transform column-wise
  results <- lapply(seq_len(n_series), function(ii) cep.mtm(x[, ii], tau))

  # Each lpa$lspec is length N (time length); rbind gives n_series x N
  lspec_mat <- do.call(rbind, lapply(results, `[[`, "lspec"))
  cep_mat   <- do.call(rbind, lapply(results, `[[`, "cep"))

  N <- ncol(lspec_mat)   # time-length
  D.hat <- as.data.frame(cep_mat, stringsAsFactors = FALSE)
  colnames(D.hat) <- paste0("C", 0:(N - 1L))

  D.hat$y <- y
  D.hat
}

#' Choose optimal number of cepstral coefficients by leave-one-out LDA CV
#'
#' @param data data.frame returned by cep.get (with columns C0..C_{N-1} and y)
#' @param mcep integer max cepstral order to evaluate
#' @return integer optimal L (index of first max CV accuracy)
Lopt.get <- function(data, mcep) {
  if (!is.data.frame(data)) stop("data must be a data.frame produced by cep.get()")
  if (!("y" %in% colnames(data))) stop("data must contain column 'y'")

  mcep <- as.integer(mcep)
  if (mcep < 1L) stop("mcep must be >= 1")

  cvK <- numeric(length = mcep)
  for (k in seq_len(mcep)) {
    # Use columns C0..Ck (which correspond to 1:(k+1) columns in data)
    cols <- colnames(data)[1:(k + 1L)]
    form <- stats::as.formula(paste("y ~", paste(cols, collapse = " + ")))
    lda_res <- MASS::lda(form, data = data, CV = TRUE)
    cvK[k] <- mean(lda_res$class == data$y)
  }
  # choose the smallest k achieving max accuracy
  Lopt <- which(cvK == max(cvK))[1L]
  Lopt
}

#' Cepstral LDA / spectral discriminant analysis
#'
#' Performs LDA on cepstral coefficients computed from M-quantile periodograms.
#' If \code{L = FALSE} (default) the function selects L by cross-validation up to \code{mcep}.
#' If \code{xNew} is provided, predicts class labels for the new set of series.
#'
#' @param y Class labels vector length = number of columns of x.
#' @param x Numeric matrix: rows = time, cols = series.
#' @param tau Quantile level.
#' @param xNew Optional numeric matrix of new series (same format as x) for prediction.
#' @param L If numeric, uses this number of cepstral coefficients; if FALSE, uses CV to select L.
#' @param mcep Maximum number of cepstral coefficients to consider when CV is requested.
#' @param cv Logical: if TRUE, compute leave-one-out CV predictions for training data.
#' @param tol Numeric tolerance passed to lda().
#' @return An object of class "ceplda" containing lda object, cep.data, Lopt, lspec (discriminant recon), and predictions.
#' @examples
#' \dontrun{
#'   res <- cep.lda(y, X, 0.5, mcep = 8)
#' }
cep.lda <- function(y, x, tau, xNew = NULL, L = FALSE, mcep = 10, cv = FALSE, tol = 1e-4) {
  # Input checks
  if (!is.matrix(x)) stop("x must be a matrix (time x series)")
  if (ncol(x) < 2L) stop("At least two series (columns) are required in x")
  if (length(y) != ncol(x)) stop("length(y) must equal number of columns of x")
  if (!is.null(xNew) && cv) stop("Cannot provide xNew and request cross-validation at the same time")

  # Build cepstral dataset
  D.hat0 <- cep.get(y, x, tau)

  # Choose L
  if (identical(L, FALSE)) {
    Lopt <- Lopt.get(D.hat0, mcep = mcep)
  } else {
    Lopt <- as.integer(L)
    if (Lopt < 1L) stop("L must be a positive integer or FALSE")
  }

  # Build formula with C0..C_{Lopt}
  cols_use <- colnames(D.hat0)[1:(Lopt + 1L)]
  form <- stats::as.formula(paste("y ~", paste(cols_use, collapse = " + ")))

  if (isTRUE(cv)) {
    C.lda <- MASS::lda(form, data = D.hat0, CV = TRUE, tol = tol)
    pre <- list(class = C.lda$class, posterior = C.lda$posterior)
    out <- list(C.lda = C.lda, cep.data = D.hat0, Lopt = Lopt, lspec = NULL, predict = pre)
  } else {
    C.lda <- MASS::lda(form, data = D.hat0, CV = FALSE, tol = tol)

    # reconstruct discriminant functions over frequency grid
    Q <- min(Lopt, length(unique(y)) - 1L)
    freq <- seq(from = 0, to = pi, by = 1 / (ncol(D.hat0) - 1L))
    dsc <- matrix(NA_real_, nrow = length(freq), ncol = Q)
    if (Q >= 1L) {
      for (q_i in seq_len(Q)) {
        dsc[, q_i] <- recon(C.lda$scaling[, q_i], freq)
      }
    }
    lspec <- list(dsc = dsc, freq = freq)

    # prediction for training or new data
    if (!is.null(xNew)) {
      if (!is.matrix(xNew)) stop("xNew must be a matrix with same layout as x")
      # build cepstral features for new data; provide dummy y of correct length
      y_dummy <- rep(NA, ncol(xNew))
      newCep <- cep.get(y_dummy, xNew, tau)
      pre <- stats::predict(C.lda, newCep, prior = C.lda$prior)
    } else {
      pre <- stats::predict(C.lda, D.hat0, prior = C.lda$prior)
    }

    out <- list(C.lda = C.lda, cep.data = D.hat0, Lopt = Lopt, lspec = lspec, predict = pre)
  }

  class(out) <- "ceplda"
  out
}

#' Reconstruct discriminant function across frequencies
#'
#' Given weight vector \code{wts} (length nw) and frequency grid \code{fqs},
#' construct Fourier basis (constant + cosines) and compute FBmat %*% wts.
#'
#' @param wts Numeric vector of weights (length nw).
#' @param fqs Numeric vector of frequencies in [0, 0.5].
#' @return Numeric vector of reconstructed discriminant values at each frequency.
recon <- function(wts, fqs) {
  nw <- length(wts)
  FBmat <- matrix(0, nrow = length(fqs), ncol = nw)
  FBmat[, 1] <- rep(1, length(fqs))
  if (nw >= 2L) {
    for (j in 2:nw) {
      FBmat[, j] <- cos(2 * pi * (j - 1L) * fqs)
    }
  }
  as.numeric(FBmat %*% wts)
}

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