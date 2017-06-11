#' simulate
#' \code{}
#' Simulation given (g,h,x10,ξ)
#' @param g function g
#' @param h function h
#' @param x0 steady state
#' @param t integer, simulation length
#' @param e vector or matrix, each row e[k, ] corresponds to ξ_k+1
#'
#' @return out: matrix of simulation output
#'
#' @export
simulate <- function(g, h, x0, t, e){
  if (is.null(t)){
  n1 <- length(x0)
  n2 <- length(g(x0))

  pre <- 1:n1
  npr <- (n1 + 1):(n1 + n2)

  out <- matrix(0, t, n1 + n2)

  out[1, pre] <- x0
  out[1, npr] <- g(x0)

  for (i in 1:(t - 1)) {
    out[i + 1, pre] <- h(out[i, pre]) + e[i, ]
    out[i + 1, npr] <- g(out[i + 1, pre])
  }
  out
  }
  else{

  }
}
