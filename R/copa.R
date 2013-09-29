## translate a matrix onto the COPA scale, i.e. 
##
## X_c_i_j <- ( X_i_j - median(X_j) ) / mad( X_j )
## 
copa <- function(X) { 
  require(matrixStats)
  X_med <- rowMedians(X, na.rm=T)
  X_mad <- rowMads(X, centers=X_med, na.rm=T)
  sweep(sweep(X, 1, X_med), 1, X_mad, '/')
}
