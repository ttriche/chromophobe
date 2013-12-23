## for ChAMP, mostly
fixBetaNAs <- function(x) {
	if(is(x, 'GenomicRatioSet')) x <- getBeta(x)
	if(any(is.na(x))) {
		require(impute)
		x <- inv.logit(impute.knn(logit(x))$data)
	} 
	return(x)
}
