\name{SIS_OAL}
\alias{SIS_OAL}
\title{Sure Independence Screening with Outcome-Adaptive Lasso (SIS-OAL)}
\description{
Performs covariate selection and estimation of the Average Treatment Effect (ATE) in high-dimensional data using a two-step procedure: Sure Independence Screening (SIS) based on conditional correlation and Outcome-Adaptive Lasso (OAL) for propensity score estimation. This method selects variables that are predictive of the outcome and adjusts for confounding when estimating treatment effects.
}
\usage{
SIS_OAL(X, A, Y, d.n = 30)
}
\arguments{
  \item{X}{A matrix or data frame of covariates (of dimension n x p).}
  \item{A}{A binary treatment indicator vector (0 or 1) of length n.}
  \item{Y}{A continuous or binary outcome variable of length n.}
  \item{d.n}{Number of covariates to retain after the SIS step (default is 30).}
}
\value{
  \item{OAL}{Estimated average treatment effect (ATE) corresponding to the lambda value that minimizes the weighted absolute mean difference (wAMD).}
  \item{mOAL}{A binary vector indicating selected variables in the propensity score model (1 if selected, 0 otherwise).}
}

\details{
This function first ranks covariates based on conditional correlation (using the \code{Causal.cor} function), then applies the Outcome-Adaptive Lasso to estimate the propensity score using a set of lambda values. The lambda that minimizes the weighted Absolute Mean Difference (wAMD) is selected. The corresponding ATE is returned, along with an indicator vector of selected covariates.
}
\references{
Shortreed, S. M., & Ertefaie, A. (2017). Outcome-adaptive lasso: Variable selection for causal inference. \emph{Biometrics}, 73(4), 1111-1122.
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
n = naug = 300

# total number of predictors
p = 1000

# threshold for screening (taken from Fan and Lv (2008))
#d.n=floor(n/log(n))
#threshold = min(d.n,p)

# set information for data augmentation
#paug=threshold

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
res=SIS_OAL(X,A,Y,d.n=30)
}
