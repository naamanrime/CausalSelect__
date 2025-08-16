
# function data gen

Data_G <- function(X,alpha_v,beta_v,bA, sig_x, linearY=TRUE,pC,pP,pI){
  Data = NULL
  X=X
  n=nrow(X)
  p=ncol(X)

  bA=bA
  pC=pC # pC: associate with both
  pP=pP # pP: associate with outcome
  pI=pI # pI: associate with treatment


  pS = p - (pC+pI+pP) # pS: associate with neither

  var.list = c( paste( "Xc", 1:pC, sep ="") , paste( "Xp", 1:pP, sep ="") ,
                 paste( "Xi", 1:pI, sep =""), paste( "XS", 1:pS, sep =""))


  #generate data

  Data = as.data.frame(X)
  names(Data) = var.list

  #simulate A treatment
  gA_x = rowSums(Data[,var.list]*matrix(alpha_v,nrow=n,ncol=length(var.list),byrow=TRUE))
  pA = expit( gA_x )
  Data$A = as.numeric( runif(n=length(pA)) < pA) # simulate A


  #simulate Y outcome
  gY_xA = rowSums(Data[,var.list]*matrix(beta_v,nrow=n,ncol=length(var.list),byrow=TRUE))
  Data$Y = gY_xA + rnorm(n=n,sd=sig_x)
  Data$Y = Data$Y + Data$A*bA

  X=as.matrix(Data[,var.list])
  A=Data$A
  Y=Data$Y
  return( list(X=X,A=A,Y=Y) )
}

