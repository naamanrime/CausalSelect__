\name{CBS}
\alias{CBS}
\title{
Ultra-high dimensional variable selection for doubly robust causal inference
}
\description{
Performs variable screening using the Causal Ball Screening (CBS) procedure and estimates the average treatment effect (ATE) using a double robust estimator.
}
\usage{
CBS(X, A, Y, d.n = 30, alpha = 0.05)
}
\arguments{
  \item{X}{A numeric matrix of covariates (n rows by p columns).}
  \item{A}{A binary treatment indicator vector of length n (0 or 1).}
  \item{Y}{A numeric outcome vector of length n.}
  \item{d.n}{Number of variables to retain after the CBS screening step (default is 30).}
  \item{alpha}{Significance level for confidence interval computation (default is 0.05).}
}
\details{
This function implements the full CBS method for confounder selection based on the ball covariance metric, followed by a double robust estimation procedure.
}
\value{
A numeric vector containing:
  \item{CBS}{Estimated average treatment effect (ATE).}
  \item{mCBS}{A binary vector indicating which variables (among the top \code{d.n}) were selected in the final model.}
}
\references{
Tang, D., Kong, D., Pan, W., Wang, L. (2023) Ultra-high dimensional variable selection for doubly robust causal inference. Biometrics, 79, 903–914. https://doi.org/10.1111/biom.13625
}
\author{
Rime Naaman
}
\note{
Requires the functions \code{Causal.cor}, \code{create_weights}, \code{wAMD_function}, and \code{DR} to be defined in the working environment. Also depends on the \pkg{glmnet} package.
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
res2=CBS(X, A, Y, d.n = 30, alpha = 0.05)
}
\keyword{causal}
\keyword{variable selection}
\keyword{double robust}
