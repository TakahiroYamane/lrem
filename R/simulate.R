#' simulate
#'
#' @param g function g
#' @param h function h
#' @param x0 steady state
#' @param t integer, simulation length
#' @param e vector or matrix, each row e[k, ]
#'
#' @return g and h functions
#'
#' @export
simulate <- function(g, h, x0, t, e){
  n1 <- length(x0)
  n2 <- length(g(x0))

  pre <- 1:n1
  npr <- (n1 + 1):(n1 + n2)

  out <- matrix(0, t, n1 + n2)

  Zero <- as.numeric(length(x0), nrow = (nrow(x0)), ncol = (ncol(x0)))
  i <- diag(nrow(e))

  out[1, pre] <- x0
  out[1, npr] <- g(x0)

  for (i in 1:(t - 1)) {
    out[i + 1, pre] <- h(out[i, pre]) + matrix(c(i, Zero)) %*% e[i, ]
    out[i + 1, npr] <- g(out[i + 1, pre])
  }
  out
}
