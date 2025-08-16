\name{CausalSelect-package}
\alias{CausalSelect-package}
\docType{package}
\title{
CausalSelect: Variable Selection for High-Dimensional Causal Inference
}

\description{
The \code{CausalSelect} package provides variable selection tools tailored for high-dimensional causal inference problems. It implements OAL, GOAL, CBS, SIS_OAL, and SIS_GOAL.
}

\details{
In high-dimensional causal inference settings where the number of covariates \eqn{p} is much larger than the sample size \eqn{n} (\eqn{p \gg n}), variable selection is crucial to control for confounding.

This package includes the following methods:
\itemize{
  \item \code{\link{OAL}}: The original Outcome-Adaptive Lasso as proposed by Shortreed and Ertefaie (2017). It uses outcome-informed variable selection in the propensity score model to estimate causal effects.
  \item \code{\link{GOAL}}: The Generalized Outcome-Adaptive Lasso, proposed by Baldé et al. (2023), improves OAL.
  \item \code{\link{SIS_GOAL}}: A screening-based extension of GOAL for ultra-high dimensional settings, using Sure Independence Screening (SIS) to preselect variables before applying GOAL.
  \item \code{\link{SIS_OAL}}: A screening-based extension of OAL that improves computational efficiency in high-dimensional scenarios.
  \item \code{\link{CBS}}: The Causal Ball Screening method, which combines conditional correlation-based screening, weighted adaptive lasso, and doubly robust estimation for treatment effect estimation.
}

The SIS-based approaches are especially useful when \eqn{p \gg n}, where matrix operations (e.g., Cholesky decomposition) become unstable or computationally prohibitive.
}

\author{
Rime Naaman and Ismaila Balde
}

\references{
Shortreed, S. M., & Ertefaie, A. (2017). \emph{Outcome-adaptive Lasso: Variable selection for causal inference}. \strong{Biometrics}, 73(4), 1111–1122.

Baldé, I., Lefebvre, G., & Ernst, O. (2023). \emph{Generalized Outcome-Adaptive Lasso for Causal Inference in High Dimensions}. arXiv preprint arXiv:2303.08882.
}

\seealso{
\code{\link{OAL}}, \code{\link{GOAL}}, \code{\link{SIS_GOAL}}, \code{\link{SIS_OAL}}, \code{\link{CBS}}
}

\keyword{package}



