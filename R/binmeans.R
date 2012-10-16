## from Martin Morgan, via R-devel, 10/03/2012
## 
binmeans <- function(x, bx) {
  if(!exists('binmean')) {
    require(inline)
    binmean <- cfunction(signature(x="numeric", bx="numeric"),
"   int nx = Rf_length(x), nb = Rf_length(bx), i, j, n;
    SEXP ans = PROTECT(NEW_NUMERIC(nb));
    double sum, *xp = REAL(x), *bxp = REAL(bx), *ansp = REAL(ans);
    sum = j = n = 0;
    for (i = 0; i < nx; ++i) {
        while (xp[i] >= bxp[j]) {
             ansp[j++] = n > 0 ? sum / n : 0;
             sum = n = 0;
        }
        n += 1;
        sum += xp[i];
    }
    ansp[j] = n > 0 ? sum / n : 0;
    UNPROTECT(1);
    return ans;
")
  }
  binmean(x, bx)
}

## test case:
##
## nx <- 4e7
## nb <- 1e3
## x <- sort(runif(nx))
## bx <- do.call(seq, c(as.list(range(x)), length.out=nb))
##
## > bx1 <- c(bx[-1], bx[nb] + 1)
## > system.time(res <- binmean(x, bx1))
##   user    system  elapsed
##   0.052    0.000    0.050
