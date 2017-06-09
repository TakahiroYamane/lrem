#' impluse
#'
#' @param g function g
#' @param h function h
#' @param x0 steady state
#' @param t integer, simulation length
#' @param e1 i-th element e[i] corresponds to a shock to i-th endogenous variable
#'
#' @return g and h functions
#'
#' @export
impluse <- function(g, h, x0, t, e1){
  n1 <- length(x0)
  n2 <- length(g(x0))

  pre <- 1:n1
  npr <- (n1 + 1):(n1 + n2)

  out <- matrix(0, t, n1 + n2)

  Zero <- as.numeric(length(x0), nrow = (nrow(x0)), ncol = (ncol(x0)))
  i <- diag(nrow(e1))

  out[1, pre] <- x0 + matrix(c(i, Zero)) %*% e1
  out[1, npr] <- g(x0)

  for (i in 1:(t - 1)) {
    out[i + 1, pre] <- h(out[i, pre])
    out[i + 1, npr] <- g(out[i + 1, pre])
  }
  out
}

