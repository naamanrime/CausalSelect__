OAL=function(X,A,Y){

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
  var.list = c(paste("X",1:p,sep=""))
  Data=NULL
  Data=as.data.frame(X)
  true_var_names<-colnames(Data)
  colnames(Data)=var.list
  Data$A=A
  Data$Y=Y


  # Normlize coviarates to have mean 0 and standard deviation 1
  temp.mean = colMeans(Data[,var.list])
  Temp.mean = matrix(temp.mean,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[,var.list] = Data[,var.list] - Temp.mean
  temp.sd = apply(Data[var.list],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[var.list] = Data[,var.list] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))


  # estimate outcome model
  y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
  lm.Y = lm(y.form,data=Data)
  beta.xy=coef(lm.Y)
  betaXY.ols = coef(lm.Y)[var.list]

  ## Want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  rownames(coeff_XA) = var.list

  ######################################################################################
  #####  Run outcome adaptive lasso for each lambda value
  ######################################################################################
  # weight model with all possible covariates included, this is passed into lasso function
  w.full.form = formula(paste("A~",paste(var.list,collapse="+")))
  for( lil in names(lambda_vec) ){
    il = lambda_vec[lil]
    ig = gamma_vals[lil]

    ### create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
    oal_pen = adaptive.lasso(lambda=n^(il),al.weights = abs(betaXY.ols)^(-ig) )
    ### run outcome-adaptive lasso model with appropriate penalty
    #logit_oal = lqa::lqa( w.full.form, data=Data, penalty=oal_pen, family=binomial(logit) )
    X=as.matrix(Data[var.list]);y=as.vector(Data$A);
    logit_oal = lqa.default( x=X, y=y, penalty=oal_pen, family=binomial(logit))
    # generate propensity score for ATE
    Data[,paste("f.pA",lil,sep="")]=expit(as.matrix(cbind(rep(1,n),Data[var.list]))%*%as.matrix(coef(logit_oal)))
    # save propensity score coefficients
    coeff_XA[var.list,lil] = coef(logit_oal)[var.list]
    # create inverse probability of treatment weights
    Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
    # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
    wAMD_vec[lil] = wAMD_function(DataM=Data,varlist=var.list,trt.var="A",
                                  wgt=paste("w",lil,sep=""),beta=betaXY.ols)$wAMD
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
  coeff_XA[,tt]

  # save the ATE corresponding the optimal lambda2 value
  OAL=ATE[tt][[1]]
  # save the propensity score coefficients corresponding the optimal lambda2 value
  mOAL=ifelse(abs(coeff_XA[,tt])> 10^(-8),1,0)
  # save wAMD value corresponding the optimal lambda2 value
  #WOAL_l[s]=wAMD_vec[tt][[1]]
  cat("Selected variables are:\n")
  print(true_var_names[which(mOAL==1)])
  return(c(OAL, mOAL))
}
