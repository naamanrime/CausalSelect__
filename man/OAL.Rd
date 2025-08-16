\name{OAL}
\alias{OAL}
\title{
Outcome-Adaptive Lasso (OAL) for High-Dimensional Causal Inference
}
\description{
Estimates the average treatment effect (ATE) using the Outcome-Adaptive Lasso (OAL)method, which performs variable selection by weighting the penalty based on the importance of covariates in the outcome model.
}
\usage{
OAL(X,A,Y)
}
\arguments{
  \item{X}{A numeric matrix or data frame of covariates with \eqn{n} rows (observations) and \eqn{p} columns (variables).}
  \item{A}{A binary treatment vector of length \eqn{n}.}
  \item{Y}{A numeric outcome vector of length \eqn{n}.}
}
\details{
This function :
\itemize{
  \item fits a linear model for the outcome to compute adaptive weights
  \item applies an adaptive lasso penalty to estimate the propensity score model
  \item computes inverse probability of treatment weights (IPTW) and ATE estimates across a grid of penalty values (\eqn{\lambda})
  \item selects the optimal \eqn{\lambda} by minimizing the weighted absolute mean difference (wAMD) across covariates
}
}
\value{
A numeric vector containing:
  \item{OAL}{Estimated average treatment effect (ATE) for the optimal lambda.}
  \item{mOAL}{A binary vector indicating selected variables (1 if selected, 0 otherwise).}
}
\references{
Shortreed, S. M., & Ertefaie, A. (2017). Outcome-adaptive lasso: Variable selection for causal inference.
  \emph{Biometrics}, 73(4), 1111–1122.
}
\author{
Rime Naaman
}
\seealso{
\code{\link{GOAL}}
}
\examples{
## Generate a multivariate normal X matrix
# set information for simulating coviariates
mean_x = 0
sig_x = 1

# pairwise correlation between covariates
rho = 0

# set number of monte carlo (MC) simulation
S=5

# sample size
n = naug = 30

# total number of predictors
p = 10

# note:  pC, pP and pI are number of confounders,
#pure predictors of outcome and pure predictors of exposure, respectively
pC = pP = pI  = 2

# pS number of spurious covariates
pS = p - (pC+pP+pI)

# list of all p variables
var.list = c(paste("Xc",1:pC,sep=""),
             paste("Xp",1:pP,sep=""),
             paste("Xi",1:pI,sep=""),
             paste("Xs",1:pS,sep=""))

# list of threshold variables
#var.list_Ball = c(paste("X",1:threshold,sep=""))

# set strength of relationship between covariates and outcome
beta_v =  c( 0.6, 0.6, 0.6, 0.6, 0, 0, rep(0,p-6) )

# Set strength of relationship between covariates and treatment
alpha_v = c( 1, 1, 0, 0, 1, 1,  rep(0,p-6) )
names(beta_v) = names(alpha_v) = var.list

# set true average treatment effect (taken from Shortreed and Ertefaie (2017))
bA = 0


Data=NULL
### simulate data
Sigma_x = matrix(rho*sig_x^2,nrow=length(var.list),ncol=length(var.list))
diag(Sigma_x) = sig_x^2
Mean_x = rep(mean_x,length(var.list))
Data = as.data.frame(MASS::mvrnorm(n = n,mu=Mean_x,Sigma = Sigma_x,empirical = FALSE))
names(Data) = var.list
X=Data



## Data generation setting
## alpha: Xc's scale is 0.2 0.2 and Xi's scale is 0.3 0.3
## so this refers that there is 2 Xc and Xi
## beta: Xc's scale is 2 2 and Xp's scale is 2 2
## so this refers that there is 2 Xc and Xp
## rest with following setup
Data_fun <- Data_G(X, alpha_v = c( 1, 1, 0, 0, 1, 1,  rep(0,p-6) )
, beta_v = c( 0.6, 0.6, 0.6, 0.6, 0, 0, rep(0,p-6) )
, bA = 0, sig_x=sig_x, linearY=TRUE,pC=2,pP=2,pI=2)

X=Data_fun$X
A=Data_fun$A
Y=Data_fun$Y
res=OAL(X,A,Y)
}
