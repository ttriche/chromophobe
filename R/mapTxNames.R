mapTxNames <- function(ids, db=NULL, idtype='TXNAME', what=c('TXID','SYMBOL')) {

  warning('FIXME: this needs MUCH better testing!')
  if(is.null(db)) {
    require(Homo.sapiens)
    db <- Homo.sapiens
  }
  ok <- !is.na(ids)
  sid <- ids[ok]
  symbolMappings <- select(db, cols=c(idtype, what), keys=sid, keytype=idtype)
  symbolMappings <- symbolMappings[ match(ids, symbolMappings[[idtype]]), ]
  hasSymbol <- which(!is.na(symbolMappings[ ids, 'SYMBOL' ]))
  syms <- symbolMappings[ ids, 'SYMBOL' ]
  syms[ -hasSymbol ] <- symbolMappings[ ids[ -hasSymbol ], what[[1]] ]
  syms

}

mapUCSCKG <- mapTxNames ## per Marc
