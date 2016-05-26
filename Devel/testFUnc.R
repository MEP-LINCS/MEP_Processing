#' Loess normalize an array using the biological residuals
#' @export
loessNorm <- function(Value,Residual,ArrayRow,ArrayColumn){
  dt <-data.table(Value=Value,Residual=Residual,ArrayRow=ArrayRow,ArrayColumn=ArrayColumn)
  lm <- loess(Residual~ArrayRow+ArrayColumn, dt, span=.7)
  dt$ResidualLoess<-predict(lm)
  dt <- dt[,ValueLoess := Value-ResidualLoess]
  return(ValueLoess = dt$ValueLoess)
}

loessNormArrays <- function(dt){
  #Loess normalize values within an array
  #Get the spot mean over all arrays
  dt <- dt[,mel :=mean(Value),by="MEP"]
  #Get the residuals from the spot mean
  dt <- dt[,Residual := Value-mel]
  #Subtract the loess model of each array's residuals from the signal
  dt <- dt[, ValueLoess:= loessNorm(Value,Residual,ArrayRow,ArrayColumn), by="BWL"]
  dt$ValueLoess[dt$k==0] <- dt$Value
  return(dt)
}

RUV3
tmp <- dt[,lapply(.SD, tmpFunct,.(Barcode)),.SDcols=signalNames]
