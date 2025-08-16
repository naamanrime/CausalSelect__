SIS_GOAL<- function(X,A,Y,d.n=30){

  # set vector of possible lambda's to try (taken from Shortreed and Ertefaie (2017))
  lambda_vec = c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
  names(lambda_vec) = as.character(lambda_vec)
  # lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  gamma_convergence_factor = 2
  # get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
  names(gamma_vals) = names(lambda_vec)

  p <- ncol(X)
  n <- nrow(X)
  d.n=floor(n/log(n))
  threshold = min(d.n,p)
  # list of threshold variables
  var.list_Ball = c(paste("X",1:threshold,sep=""))

  #set threshold for screening.
  ballcor<-rep(NA, p)
  for (j in 1:p){
    # calculate conditional ball covariance for each variable.
    ballcor[j]<-Causal.cor(X[,j],Y,A)
  }

  # screening procedure
  ballorder<-order(ballcor, decreasing=TRUE)
  # select the top threshold number of variables
  ballsisselectedindex<-ballorder[1:threshold]
  ballsisselectedindex = ballsisselectedindex[order(ballsisselectedindex)]
  weight = ballcor[ballsisselectedindex]

  # the data matrix after screening
  Data=NULL
  Data = X[,ballsisselectedindex]
  true_var_names<-colnames(Data)
  Data = as.data.frame(Data)
  names(Data) = var.list_Ball
  Data$A = A
  Data$Y = Y

  # normlize coviarates to have mean 0 and standard deviation 1
  temp.mean = colMeans(Data[,var.list_Ball])
  Temp.mean = matrix(temp.mean,ncol=length(var.list_Ball),nrow=nrow(Data),byrow=TRUE)
  Data[,var.list_Ball] = Data[,var.list_Ball] - Temp.mean
  temp.sd = apply(Data[var.list_Ball],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list_Ball),nrow=nrow(Data),byrow=TRUE)
  Data[var.list_Ball] = Data[,var.list_Ball] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))

  # estimate outcome model
  y.form = formula(paste("Y~A+",paste(var.list_Ball,collapse="+")))
  lm.Y = lm(y.form,data=Data)
  betaXY = coef(lm.Y)[var.list_Ball]
  summary(lm.Y)
  betaXY

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%                 create Augmented data      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  # set the possible lambda2 value (taken from Zou and Hastie (2005))
  S_lam=c(0,10^c(-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1))
  WM_N=WM_P=M_N=M_P=S_wamd=rep(NA,length(S_lam))
  M_mat=matrix(NA,length(S_lam),threshold)
  for (m in 1:length(S_lam)) {
    # augmented A and X
    lambda2=S_lam[m]
    q=threshold
    I=diag(1,q,q)
    Iq=sqrt(lambda2)*I
    Anq=c(Data$A,rep(0,q))
    Xnq=matrix(NA,n+q,q)
    X_th=Data[,var.list_Ball]
    for (j in 1:q){
      Xnq[,j]=c(X_th[,j],Iq[,j])
    }
    newData=as.data.frame(Xnq)
    names(newData)=var.list_Ball
    newData$A=Anq
    n.q=n+q

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%%%%%%%%%%%%%% run GOAL based on PIRLS                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # want to save ATE, wAMD and propensity score coefficients for each lambda value
    ATE = wAMD_vec = rep(NA, length(lambda_vec))
    names(ATE) = names(wAMD_vec) = names(lambda_vec)
    coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list_Ball),ncol=length(lambda_vec)))
    names(coeff_XA) = names(lambda_vec)
    rownames(coeff_XA) = var.list_Ball

    # weight model with all possible covariates included, this is passed into lasso function
    for( lil in names(lambda_vec) ){
      il = lambda_vec[lil]
      ig = gamma_vals[lil]

      # create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
      oal_pen = adaptive.lasso(lambda=n.q^(il),al.weights = abs(betaXY)^(-ig) )
      # run outcome-adaptive lasso model with appropriate penalty
      X=as.matrix(newData[var.list_Ball]);y=as.vector(newData$A);
      lq22=lqa.def(x=X, y=y, family = binomial,penalty = oal_pen)
      # generate propensity score for ATE
      Data[,paste("f.pA",lil,sep="")]=expit(as.matrix(cbind(rep(1,n),Data[var.list_Ball]))%*%as.matrix((1+lambda2)*coef(lq22)))
      # create inverse probability of treatment weights for ATE
      Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
      #save propensity score coef
      coeff_XA[var.list_Ball,lil] = (1+lambda2)*coef(lq22)[var.list_Ball]
      # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
      wAMD_vec[lil] = wAMD_function(DataM=Data,varlist=var.list_Ball,trt.var="A",
                                    wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD

      # save ATE estimate for this lambda value
      ATE[lil] = ATE_est(fY=Data$Y,fw=Data[,paste("w",lil,sep="")],fA=Data$A)
    } # close loop through lambda values

    # print out wAMD for all the lambda values tried
    wAMD_vec
    # find the lambda value that creates the smallest wAMD
    tt = which.min(wAMD_vec)
    # print out ATE corresponding to smallest wAMD value
    ATE[tt]
    # print out the coefficients for the propensity score that corresponds with smalles wAMD value
    GOAL.PIRLS=coeff_XA[,tt]
    GOAL.PIRLS.ate=ATE[tt][[1]]

    M_P[m]=GOAL.PIRLS.ate
    WM_P[m]=wAMD_vec[tt][[1]]
    M_mat[m,]=ifelse(abs(coeff_XA[,tt])> 10^(-8),1,0)

  } # end for for lambda2
  # final ATE and coe
  ptt= which.min(WM_P)
  GOAL=M_P[ptt]
  mGOAL=ifelse(abs(M_mat[ptt,])> 10^(-8),1,0)
  cat("Selected variables are:\n")
  print(true_var_names[which(mGOAL==1)])
  return(c(GOAL,mGOAL))
}
