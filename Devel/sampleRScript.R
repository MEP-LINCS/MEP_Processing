# R in C Sample Script
#Mark Dane 1/2016

transformArrays <- function(x,y, fileName="RPlot.pdf"){
  xp <- x*5
  yp <- y*10
  pdf(file=fileName)
  plot(x=x,y=xp)
  dev.off()
  
  return(list(xPrime=xp, yPrime=yp))
}

fList <- transformArrays(x=1:5, y=10:20, fileName="Sample.pdf")
