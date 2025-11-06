#' Discretized beta distribution
#' R implementation of the RevBayes function fnDiscretizeBeta
#' @source https://gist.githubusercontent.com/ms609/883632d10d4d80ea5391cee9c47071fc/raw/505b31f322060dd2cdde8010ae74708751ee862e/fnDiscretizeBeta.R
#' Martin R. Smith
#' 2024-12-05
#' @export
incompleteBeta <- function(a, b, x) {
  tol <- sqrt(.Machine[["double.eps"]])
  
  if (x <= 0) {
    0
  } else if (1 <= x) {
    1
  } else {
    # Change tail if necessary and determine S
    psq <- a + b
    if (a < psq * x) {
      xx <- 1 - x
      cx <- x
      pp <- b
      qq <- a
      indx <- TRUE
    } else {
      xx <- x
      cx <- 1 - x
      pp <- a
      qq <- b
      indx <- FALSE
    }
    
    term <- 1
    i <- 1L
    value <- 1
    ns <- floor(qq + cx * psq)
    
    # use Soper's reduction formulae
    rx <- xx / cx
    
    temp <- qq - i
    if (ns == 0) {
      rx <- xx
    }
    
    for (it in 1:1001) {
      if (it == 1001) return(Inf)
      
      term <- term * temp * rx / (pp + i)
      value <- value + term
      temp <- abs(term)
      if (temp <= tol && temp <= tol * value) {
        break
      }
      i <- i + 1
      ns <- ns - 1
      if (0 <= ns) {
        temp <- qq - i
        if (ns == 0) {
          rx <- xx
        }
      } else {
        temp <- psq
        psq <- psq + 1.0
      }
    }
    
    # Finish calculation
    value <- value * exp(pp * log(xx) + (qq - 1.0) * log(cx) - lbeta(a, b) - log(pp))
    
    # Return:
    if (indx) {
      1.0 - value
    } else {
      value
    }
  }
}

#' @rdname incompleteBeta
#' @export
fnDiscretizeBeta <- function(alpha, beta, nCats, median = FALSE) {
  factor <- alpha / sum(alpha, beta) * nCats
  if (median) {
    q <- qbeta(1:nCats / nCats - (1/(nCats * 2)), alpha, beta)
    factor * q / sum(q)
  } else {
    q <- qbeta(1:(nCats - 1) / nCats, alpha, beta)
    q2 <- vapply(q, function(x) incompleteBeta(alpha + 1, beta, x), double(1))
    factor * (c(q2, 1) - c(0, q2))
  }
}
