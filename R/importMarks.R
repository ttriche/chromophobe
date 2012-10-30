import.marks <- function(file) {

  require(utils)
  w <- rep('character', 4)
  marks <- suppressWarnings(read.delim(file, header=F, skip=1, fill=T, colCl=w))
  names(marks) <- c('Cell','Mark','Treatment','Control')
  return(marks)

}
