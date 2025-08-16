\name{Data_G}
\alias{Data_G}
\title{Simulate Treatment and Outcome Data from Covariates}
\description{
This function generates simulated data for estimating the average treatment effect (ATE).
}
\usage{
Data_G(X, alpha_v, beta_v, bA, sig_x, linearY = TRUE, pC, pP, pI)
}
\arguments{
  \item{X}{A matrix of covariates (\eqn{n \times p}).}
  \item{alpha_v}{A vector of coefficients for the treatment assignment model.}
  \item{beta_v}{A vector of coefficients for the outcome model.}
  \item{bA}{True treatment effect to be added to the outcome.}
  \item{linearY}{Logical. If \code{TRUE}, generates a linear outcome. (Currently always linear).}
  \item{sig_x}{Standard deviation of the covariates in the simulated data.}
  \item{pC}{Number of covariates associated with both treatment and outcome.}
  \item{pP}{Number of covariates associated only with the outcome.}
  \item{pI}{Number of covariates associated only with the treatment.}
}
\value{
A list containing:
  \item{X}{Matrix of covariates}
  \item{A}{Binary treatment vector}
  \item{Y}{Continuous outcome vector affected by covariates and treatment effect.}
}
\author{Rime Naaman and Ismaila Baldé}
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
p = 100

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
}


