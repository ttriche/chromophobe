write.delim <- function(df, file, quote=FALSE, row.names=FALSE, sep="\t", ...) {
  # originally from the 'caroline' package on CRAN
  write.table(df, file, quote=quote, row.names=row.names, sep=sep, ...)
}
