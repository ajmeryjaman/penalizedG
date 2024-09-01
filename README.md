## Implemention of "Penalized G-estimation for effect modifier selection in a structural nested mean model for repeated outcomes"

The repository contains the following folders:

Folder | Description
--- | ---
R | Contains the source codes
man | Contains the documentation of each function used

The R folder contains the following files:

File | Description
--- | ---
[penalizedG.r](penalizedG.r) | Contains the main function penalizedG() which implements our method
[otherFUNCTIONS.r](otherFUNCTIONS.r) | Contains other required functions called in the main function

The function penalizedG() in [penalizedG.r](penalizedG.r) performs penalized G-estimation for a given longitudinal data, a specific working correlation structure and a sequence of tuning parameters, and returns result under the optimal value (selected by a double-robust information  criterion) of the tuning parameter in the given range. Currently, the function allows a continuous outcome and a binary treatment/exposure. The outcome, the exposure and the potential confounders, all can be time-varying, but the potential confounders should be continuous or binary.
#' @param data A data frame containing the variables in longitudinal format.

Please see the example given in [penalizedG.r](penalizedG.r) to generate a longitudinal data set and perform the estimation. Or, do the following:

#### R commands for installing and using our package

library(devtools) # if already installed, otherwise need to install it first

install_github("ajmeryjaman/penalizedG")

library(penalizedG)

?penalizedG
