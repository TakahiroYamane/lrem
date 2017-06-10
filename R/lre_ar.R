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
  A <- as.matrix(A)
  E <- as.matrix(E)
  B <- as.matrix(B)
  Phi <- as.matrix(Phi)

　pren <- nx + length(Phi)

  Zero1 <- matrix(0, nrow = nrow(B), ncol = ncol(B))
  Zero2 <- matrix(0, nrow = nrow(Phi), ncol = ncol(A))
  i <- diag(nrow(Phi))

  AA <- cbind(rbind(Phi, B), rbind(Zero2, A))
  EE <- cbind(rbind(i, Zero1), rbind(Zero2, A))

  ret <- QZ::qz(AA, EE)
  ord <- abs(ret$ALPHA / ret$BETA) <= 1
  ret2 <- QZ::qz.dtgsen(ret$S, ret$T, ret$Q, ret$Z, select = ord)
  Z1S <- ret2$Z[1:pren,1:pren]
  Z2S <- ret2$Z[(pren + 1):(nrow(ret2$Z))]
  SSS <- ret2$T[1:pren, 1:pren]
  TSS <- ret2$S[1:pren, 1:pren]

  g <- function(x0){
    Z2S %*% solve(Z1S) %*% x0
  }
  h <- function(x0){
    Z1S %*% solve(SSS) %*% TSS %*% solve(Z1S) %*% x0
  }
  list(g = g,h = h)
}
