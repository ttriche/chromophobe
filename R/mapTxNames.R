mapTxNames <- function(ids, db=NULL) { 

  if(is.null(db)) {
    require(Homo.sapiens)
    db <- Homo.sapiens
  }
  symbolMappings <- select(db, cols=c('TXNAME','ENTREZID','SYMBOL'), 
                           keys=symbol(txtr), keytype='TXNAME')

  hasSymbol <- !is.na(symbolMappings$SYMBOL[ ids ]) 
  ids[hasSymbol] <- symbolMappings$SYMBOL[ ids[hasSymbol] ]
  if(length(hasSymbol)<length(ids)) warn('Symbols were not found for some IDs')
  return(ids)

}

mapUCSCKG <- mapTxNames ## per Marc
