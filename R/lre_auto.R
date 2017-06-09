#' lre_auto
#'
#' @param A coeffients of previous period
#' @param E coeffients of latter period
#' @param nx number od predetermined variables
#'
#' @return g and h functions
#'
#' @export
lre_auto <- function(A, E, nx) {
  if (is.null(E)) {
    lre_auto_bk(A, nx)
  } else {
    lre_auto_klein(A, E, nx)
  }
}

lre_auto_bk <- function(A, nx) {
  npr <- nx
  A <-as.matrix(A)
  Asch <- Matrix::Schur(A)
  Asch2 <- QZ::qz.dtrsen(Asch$T, Asch$Q, abs(Asch$EValues) <= 1)
  Q1S <- Asch2$Q[1:npr,1:npr]
  Q2S <- Asch2$Q[(npr+1):(nrow(Asch2$Q)),1:npr]

  g <- function(x0){
    Q2S %*% solve(Q1S) %*% x0
  }
  h <- function(x0){
    (A[1:npr, 1:npr] %*% x0) + (A[1:npr, (npr+1):ncol(A)] %*% Q2S %*% solve(Q1S)) %*% x0
  }
  list(g = g, h = h)
}

lre_auto_klein <- function(A, E, nx) {
  npr <- nx
  A <- as.matrix(A)
  E <- as.matrix(E)
  ret <- QZ::qz(A, E)
  abs(ret$ALPHA / ret$BETA)
  ord <- abs(ret$ALPHA / ret$BETA) <= 1
  ret2 <- QZ::qz.dtgsen(ret$S, ret$T, ret$Q, ret$Z,select = ord)
  Z1S <- ret2$Z[1:npr,1:npr]
  Z2S <- ret2$Z[(npr+1):(nrow(ret2$Z))]
  SSS <- ret2$T[1:npr, 1:npr]
  TSS <- ret2$S[1:npr, 1:npr]

  g <- function(x0){
    Z2S %*% solve(Z1S) %*% x0
  }
  h <- function(x0){
    Z1S %*% solve(SSS) %*% TSS %*% solve(Z1S) %*% x0
  }
  list(g = g, h = h)
}


