\name{SIS_GOAL}
\alias{SIS_GOAL}
\title{SIS_GOAL: ATE Estimation with SIS and Penalized Logistic Regression}
\description{
This function implements the SIS + GOAL procedure, combining Sure Independence Screening (SIS) using causal ball correlation and outcome-adaptive lasso (GOAL) with data augmentation for estimating the Average Treatment Effect (ATE) in high-dimensional settings.
}
\usage{
SIS_GOAL(X, A, Y,d.n=30)
}
\arguments{
  \item{X}{covariate matrix n x p.}
  \item{A}{Binary treatment vector  (1 for treated and 0 otherwise.}
  \item{Y}{Outcome vector}
  \item{d.n}{An optional integer specifying the screening threshold (default is 30}

}
\details{
 \code{SIS_GOAL} first applies causal ball correlation screening to select the top \code{d.n} covariates most associated
  with the outcome and treatment. Then, for a grid of regularization parameters, it estimates the propensity score using
  outcome-adaptive lasso with data augmentation (GOAL). The ATE is estimated using inverse probability of treatment weighting (IPTW).
  The final model is selected based on the configuration that minimizes the weighted absolute mean difference (wAMD).
}
\value{
A vector containing:
  \item{GOAL}{Estimated ATE from the selected model.}
  \item{mGOAL}{Binary vector indicating selected variables (1 = selected).}
}
\references{
Baldé, I., Shortreed, S. M., and Ertefaie, A. (2023). \emph{Outcome-adaptive penalized estimation for causal inference in high dimensions}. JASA.
}
\author{
Rime Naaman
}
\note{
The main difference between \code{SIS_GOAL} and \code{GOAL} is that \code{SIS_GOAL} first performs a variable screening step to select a subset of important variables before applying the GOAL method.

This screening helps handle situations where the number of variables \code{p} is much larger than the sample size \code{n} by reducing dimensionality before penalized estimation.

In contrast, \code{GOAL} applies the penalized estimation directly without preliminary screening.
}

\seealso{
  This function internally uses \code{ATE_est}, \code{Causal.cor}, \code{create_weights}, \code{wAMD_function}, \code{lqa.upd}, \code{lqa.def}, \code{expit}
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
res=SIS_GOAL(X,A,Y,d.n=30)
}
