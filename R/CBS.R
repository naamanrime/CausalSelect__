CBS<- function(X,A,Y,d.n=30,alpha=0.05){

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

  threshold = min(d.n,p)


  #threshold = min(d.n,p)
  #set threshold for screening.
  ballcor<-rep(NA, p)
  X=as.matrix(X)
  for (j in 1:p){
    # calculate conditional ball covariance for each variable.
    ballcor[j]<-Causal.cor(X[,j],Y,A)
  }


  var.list_Ball = c(paste("X",1:threshold,sep=""))
  var.list = c(paste("X",1:p,sep=""))

  # set lambda values for grid search.

  lambda = log(max(p,n))^0.75/n^0.5
  lambda_val = seq(lambda*0.1,lambda*10,length.out= 10)
  names(lambda_val) = paste0("lambda",1:10)

  # set gamma values for grid search.

  gamma_vals = seq(3,20,1)
  names(gamma_vals) = paste0("gamma",gamma_vals)

  # screening procedure
  ballorder<-order(ballcor, decreasing=TRUE)

  # select the top threshold number of variables
  ballsisselectedindex<-ballorder[1:threshold]
  ballsisselectedindex = ballsisselectedindex[order(ballsisselectedindex)]
  weight = ballcor[ballsisselectedindex]

  # the data matrix after screening
  Data = X[,ballsisselectedindex]
  true_var_names<-colnames(Data)
  Data = as.data.frame(Data)
  names(Data) = var.list_Ball
  Data$A = A
  Data$Y = Y


  # centerlize, standardlize
  temp.mean = colMeans(Data[,var.list_Ball])
  Temp.mean = matrix(temp.mean,ncol=length(var.list_Ball),nrow=nrow(Data),byrow=TRUE)
  Data[,var.list_Ball] = Data[,var.list_Ball] - Temp.mean
  temp.sd = apply(Data[var.list_Ball],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list_Ball),nrow=nrow(Data),byrow=TRUE)
  Data[var.list_Ball] = Data[,var.list_Ball] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))


  # weight for each variable for refined selection
  betaXY = weight
  betaXY = weight/max(weight)

  ## Want to save wAMD and propensity score coefficients for
  ## each lambda and gamma value

  wAMD_vec = rep(NA, length(lambda_vec)*length(gamma_vals))
  names(wAMD_vec) = paste0( rep(names(lambda_vec),each = length(gamma_vals)),names(gamma_vals))
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list_Ball),ncol=length(wAMD_vec)))
  names(coeff_XA) = names(wAMD_vec)
  rownames(coeff_XA) = var.list_Ball

  ######################################################################################
  #####  Run outcome adaptive lasso for each lambda and gamma value ####################
  ######################################################################################
  # weight model with all possible covariates included, this is passed into lasso function


  for( lil in names(lambda_val) )
    for(mim in names(gamma_vals)){
      il = lambda_val[lil]
      ig = gamma_vals[mim]
      fitx = as.matrix(Data[,1:threshold])
      alasso <- glmnet(x = fitx, y = Data$A,
                       type.measure = "class",
                       family = "binomial",
                       alpha = 1,
                       penalty.factor = c(abs(betaXY)^(-ig)),
                       lambda = il)

      # calculate propensity score
      prob = predict(alasso,newx = fitx)
      prob = exp(prob)/(1+exp(prob))
      Data[,paste("f.pA",paste0( lil,mim),sep="")] = prob
      # save propensity score coefficients
      coeff_XA[var.list_Ball,paste0( lil,mim)] = coef.glmnet(alasso,s = alasso$lambda.min)[var.list_Ball,]
      # create inverse probability of treatment weights
      Data[,paste("w",paste0( lil,mim),sep="")] = create_weights(fp=Data[,paste("f.pA",paste0( lil,mim),sep="")],fA=Data$A)
      # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
      wAMD_vec[paste0( lil,mim)] = wAMD_function(DataM=Data,varlist=var.list_Ball,trt.var="A",
                                                 wgt=paste("w",paste0( lil,mim),sep=""),beta=betaXY)$wAMD
      # save ATE estimate for this lambda value
    } # close loop through lambda values

  # find the target (gamma,lambda) with smallest wAMD score.
  tt = which.min(wAMD_vec)
  #save the estimated propensity score model
  fitted.ps = Data[,paste("f.pA",names(tt),sep="")]
  #outcome regression
  {
    # use lasso to fit the treated group.
    X1 = X[A==1,]
    Y1 = Y[A==1]
    X0 = X[A==0,]
    Y0 = Y[A==0]

    lasso<-cv.glmnet(x = X1, y = Y1,
                     type.measure = "mse",
                     ## K = 10 is the default.
                     nfolds = 10,
                     alpha = 1)
    coeflasso1 <- coef.glmnet(lasso,lasso$lambda.min)
    index1 = which(coeflasso1[-1]!=0)
    X1.fit = X1[,index1]

    data1 = data.frame(Y1,X1.fit)
    names(data1) = c("Y",paste0("v",1:length(index1)))[1:ncol(data1)]

    ormodel1 <- lm(Y~.,data1)

    data.fit = data.frame(Y,X[,index1])
    names(data.fit) = c("Y",paste0("v",1:length(index1))) [1:ncol(data1)]
    # save predict for the treated group.
    orfit1 = predict.lm(ormodel1, newdata = data.fit)

    # use lasso to fit the controlled group.
    lasso<-cv.glmnet(x = X0, y = Y0,
                     type.measure = "mse",
                     ## K = 10 is the default.
                     nfolds = 10,
                     alpha = 1)
    coeflasso0 <- coef.glmnet(lasso,lasso$lambda.min)
    index0 = which(coeflasso0[-1]!=0)
    X0.fit = X0[,index0]

    data0 = data.frame(Y0,X0.fit)
    names(data0) = c("Y",paste0("v",1:length(index0)))[1:ncol(data0)]

    ormodel0 <- lm(Y~.,data0)
    data.fit = data.frame(Y,X[,index0])
    names(data.fit) = c("Y",paste0("v",1:length(index0))) [1:ncol(data0)]
    #save the estimated data for the controlled group.
    orfit0 = predict.lm(ormodel0, newdata = data.fit)
  }
  # get the double robust estimate.
  result = DR(fitted.ps,orfit1,orfit0,A,Y)

  # get the point estimate and variance of the resulting estimator.
  result = c(mean(result),var(result))
  c("point estimate" = result[1],
    "lower bound" = result[1]-qnorm(1-alpha/2)*sqrt(result[2]/n),
    "upper bound" = result[1]+qnorm(1-alpha/2)*sqrt(result[2]/n),
    "variance" =  sqrt(result[2]/n))


  CBS=result[1]
  mCBS=ifelse(abs(coeff_XA[var.list_Ball,tt])> 10^(-8),1,0)
  cat("Selected variables are:\n")
  print(true_var_names[which(mCBS==1)])
  return(c(CBS,mCBS))

}
