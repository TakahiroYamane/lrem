#' lre_ar
#'
#' @param A coeffients of previous period
#' @param E coeffients of latter period
#' @param B matrixB
#' @param Phi matrixPhi
#' @param nx number od predetermined variables
#'
#' @return g and h functions
#'
#' @export
lre_ar <- function(A, E, B, Phi, nx){
  npr <- nx
  A <- as.matrix(A)
  E <- as.matrix(E)
  B <- as.matrix(B)
  Phi <- as.matrix(Phi)

  Zero <- matrix(0, nrow = nrow(Phi), ncol = ncol(A))
  i <- diag(nrow(Phi))
  A2 <- matrix(c(Phi, B, Zero, A), nrow = (nrow(Phi) + nrow(B)), ncol = (ncol(B) + ncol(A)))
  E2 <- matrix(c(i, Zero, Zero, A), nrow = (nrow(Phi) + nrow(B)), ncol = (ncol(B) + ncol(A)))

  ret <- QZ::qz(A2, E2)
  abs(ret$ALPHA / ret$BETA)
  ord <- abs(ret$ALPHA / ret$BETA) <= 1
  ret2 <- QZ::qz.dtgsen(ret$S, ret$T, ret$Q, ret$Z,select = ord)
  Z1S <- ret2$Z[1:npr,1:npr]
  Z2S <- ret2$Z[(npr+1):(nrow(ret2$Z))]
  SSS <- ret2$T[1:npr, 1:npr]
  TSS <- ret2$S[1:npr, 1:npr]

  g <- function(u0, x0){
    Z2S %*% solve(Z1S) %*% matrix(c(u0, x0))
  }
  h <- function(u0, x0){
    Z1S %*% solve(SSS) %*% TSS %*% solve(Z1S) %*% matrix(c(u0, x0))
  }
  list(g = g,h = h)
}
