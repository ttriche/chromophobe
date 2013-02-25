mapTxNames <- function(ids, db=NULL, idtype='TXNAME', mapping='SYMBOL') { 

  if(is.null(db)) {
    require(Homo.sapiens)
    db <- Homo.sapiens
  }
  symbolMappings <- select(db, cols=c(idtype, mapping), 
                           keys=ids, keytype=idtype)
  symbolMappings <- symbolMappings[ match(ids, symbolMappings[[idtype]]), ]
  hasSymbol <- !is.na(symbolMappings[[mapping]][ ids ]) 
  ids[hasSymbol] <- symbolMappings[[mapping]][ ids[hasSymbol] ]
  if(length(hasSymbol) < length(ids)) warn('Symbols not found for some IDs')
  return(ids)

}

mapUCSCKG <- mapTxNames ## per Marc
