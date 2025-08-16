GOAL<- function(X,A,Y){

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
  n.p=n+p
  var.list = c(paste("X",1:p,sep=""))
  Data=NULL
  Data=as.data.frame(X)
  true_var_names<-colnames(Data)
  colnames(Data)=var.list
  Data$A=A
  Data$Y=Y
  # normalize covariates to have mean 0 and standard deviation 1
  temp.mean = colMeans(Data[,var.list])
  Temp.mean = matrix(temp.mean,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[,var.list] = Data[,var.list] - Temp.mean
  temp.sd = apply(Data[var.list],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[var.list] = Data[,var.list] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
  w.full.form = formula(paste("A~",paste(var.list,collapse="+")))

  # estimate outcome model
  y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
  lm.Y = lm(y.form,data=Data)
  betaXY = coef(lm.Y)[var.list]

  ## want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  rownames(coeff_XA) = var.list

  ## set the possible lambda2 value (taken from Zou and Hastie (2005))
  S_lam=c(0,10^c(-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1))
  ## want to save ATE, wAMD and propensity score coefficients for each lambda2 value
  WM_P=M_P=S_wamd=rep(NA,length(S_lam))
  M_mat=matrix(NA,length(S_lam),p)
  for (m in 1:length(S_lam)) {
    ## create augmented A and X
    lambda2=S_lam[m]
    I=diag(1,p,p)
    Ip=sqrt(lambda2)*I
    Anp=c(Data$A,rep(0,p))
    Xnp=matrix(0,n+p,p)
    X=Data[,var.list]
    for (j in 1:p){
      Xnp[,j]=c(X[,j],Ip[,j])
    }
    newData=as.data.frame(Xnp)
    names(newData)=var.list
    newData$A=Anp

    ## want to save ATE, wAMD and propensity score coefficients for each lambda value
    ATE_try=ATE = wAMD_vec =coeff_XA=NULL;
    ATE_try=ATE = wAMD_vec = rep(NA, length(lambda_vec))
    names(ATE) = names(wAMD_vec) = names(lambda_vec)
    coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
    names(coeff_XA) = names(lambda_vec)
    rownames(coeff_XA) = var.list

    ## run GOAL based on PIRLS
    # weight model with all possible covariates included, this is passed into lasso function
    for( lil in names(lambda_vec) ){
      il = lambda_vec[lil]
      ig = gamma_vals[lil]

      # create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
      oal_pen = adaptive.lasso(lambda=n.p^(il),al.weights = abs(betaXY)^(-ig) )
      # run outcome-adaptive lasso model with appropriate penalty
      X=as.matrix(newData[var.list]);y=as.vector(newData$A);
      lq22=lqa.def(x=X, y=y, family = binomial,penalty = oal_pen)

      # generate propensity score for ATE
      Data[,paste("f.pA",lil,sep="")]=expit(as.matrix(cbind(rep(1,n),Data[var.list]))%*%as.matrix((1+lambda2)*coef(lq22)))
      # create inverse probability of treatment weights for ATE
      Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
      # save propensity score coef
      coeff_XA[var.list,lil] = (1+lambda2)*coef(lq22)[var.list]
      # estimate weighted absolute mean different over all covariates using this lambda to generate weights
      wAMD_vec[lil] = wAMD_function(DataM=Data,varlist=var.list,trt.var="A",
                                    wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD

      # save ATE estimate for this lambda value
      ATE[lil] = ATE_est(fY=Data$Y,fw=Data[,paste("w",lil,sep="")],fA=Data$A)
    } # close loop through lambda values

    # print out wAMD for all the lambda values evaluated
    wAMD_vec
    # find the lambda value that creates the smallest wAMD
    tt = which.min(wAMD_vec)
    # print out ATE corresponding to smallest wAMD value
    ATE[tt][[1]]
    # save the coefficients for the propensity score that corresponds to smallest wAMD value
    GOAL.i.c=coeff_XA[,tt]
    # check which covariates are selected
    M_mat[m,]=ifelse(abs(coeff_XA[,tt])> 10^(-8),1,0)
    # save the ATE corresponding to smallest wAMD value
    GOAL_est1=ATE[tt][[1]]
    GOAL.i.ate=GOAL_est1
    M_P[m]=GOAL.i.ate
    # save the smallest wAMD value
    WM_P[m]=wAMD_vec[tt][[1]]

  }

  ptt= which.min(WM_P)
  GOAL=M_P[ptt]
  mGOAL=ifelse(abs(M_mat[ptt,])> 10^(-8),1,0)
  cat("Selected variables are:\n")
  print(true_var_names[which(mGOAL==1)])

  return(c(GOAL,mGOAL))
}
