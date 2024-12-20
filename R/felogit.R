#' This function implements the estimators of the bounds on the AME/ATE proposed in DDL.
#'
#' The function summary_felogit applied to the output of this felogit function prints
#' the table containing the estimation results.
#'
#' @param data is an environment variable containing the data. It can either be
#' in long format or in wide format. If in long format, data should be a data
#' frame with individual indexes in the first columns, periods in the second
#' column, and in the remaining columns the value of each variable at the
#' corresponding period for the corresponding individual. All columns should
#' be named. If in wide format, data should be a list or environment variable
#' containing:
#' - data$Y is a matrix of size n x Tmax containing the values of the dependent variable Y
#' for each individual at each period. NAs represent missing observations.
#' - data$X is an array of size n x Tmax x dimX containing the values of the covariates X
#' for each individual at each period. NAs represent missing observations.
#' - data$clusterIndexes is a vector of size n x 1 that specifies the cluster each
#' observation pertains to. If it does not exist, the function enforces the default
#' setting of i.i.d. observations - the parameter takes value 1:n so that each
#' observation is its own cluster.
#' @param formul (default NULL) is a formula used when the data is in long format
#' to specify which variable to use. It must be of the form formula("Y ~ X_1 + X_2")
#' where Y is the name of the binary variable of interest and X_1, X_2, etc are
#' the variables to use in the logit model. Set formul to NULL to indicate the
#' data is already in wide format.
#' @param Option (default "quick") Estimation method to be used. If "quick" or
#' "outer" (case-insensitive) the outer bounds are computed. Otherwise, the sharp
#' bounds are computed. If any of the variable with respect to which the effect is
#' being measured is binary (so that an ATE must be computed), we always switch to
#' the quick method for efficiency reasons. In general, we recommend using the
#' outer bounds method if the number of covariates is at least three, if the
#' number of periods observed is four or more, or if the sample size is small
#' (less than 500) or large (more than 10^4).
#' @param compute_X (default "all") is a vector containing all the variables with
#' respect to which the AME/ATE must be computed. It can either contain the
#' variable names, as given in the arguments formul or varX, or their rank, e.g.
#' 3 for the third variable appearing in formul/varX or, for lack theoreof, the third
#' column in data$X.
#' @param compute_T (default "all") is a vector containing all periods at which
#' the AME/ATE must be computed. Alternatively, it can be "all", in which case the
#' AME/ATE will be computed successively at every period and, on top of that, the
#' average AME/ATE across all periods will also be computed using the function
#' compute_average_AMTE. Also note that non-positive values will be counted
#' backwards from the last period at which each individual is observed, as in
#' an event-study. If NULL, the first period is used. Periods MUST be specified
#' by their rank (1 for first period, etc) and NOT by value (e.g. 1980 for the
#' year).
#' @param cluster is clustering
#' @param alpha (default 0.05) desired asymptotic level of the estimated
#' confidence intervals
#' @param CIOption (default "CI2") When the outer bounds method is being used,
#' specifies which confidence interval should be used. If "CI2", the CI2
#' confidence interval is being used (see DDL, section 4.2), otherwise the CI3
#' confidence interval will be used (see DDL, appendix C). We recommend using
#' CI3 only if the user suspects the FE logit model may be a severely
#' misspecified model for the data.
#' @param nbCores (default 4) number of cores to be used for parallel computing,
#' to speed up the estimation of the sharp bounds.
#' @param varY (default "Y") for data already in wide format, the name of the
#' binary variable of interest in the data list/environment.
#' @param varX (default NULL) for data already in wide format, the name to use
#' for each of the variables given by the slices of the array data$X along
#' the third dimension. dimnames(data$X) along the third dimension can also be
#' used. If varX is NULL and no name is given in dimnames(data$X), X_1, X_2, etc
#' will be used.
#'
#'
#'
#' @return A list containing:
#'  - summary: a dataframe containing the estimation results,
#'  - n: the number of used individuals,
#'  - ndiscard: the number of discarded individuals,
#'  - Tmax: the total number of distinct periods of observed,
#'  - vardiscard: the label of the discarded variables,
#'  - formul: the formula used (implicitly deduced if the input was NULL),
#'  - alpha: the level used for the confidence intervals,
#'  - Option: the method used,
#'  - summary_CMLE : a dataframe containing the estimation results of the CMLE.
#'
#' @export
#'
#' @examples
#'library(pglm)
#'data("UnionWage", package = "pglm")
#'
#'UnionWage$union <- UnionWage$union == "yes"
#'UnionWage$rural <- UnionWage$rural == "yes"
#'UnionWage$black <- UnionWage$com == "black" # used as test for discarded variable because constant
#'UnionWage$NorthEast <- UnionWage$region == "NorthEast"
#'sub <- UnionWage[UnionWage$year < 1986,]
#'
#'formul <- formula("union ~ exper + married + black")
#'output <- felogit(data = sub, formul = formul, Option = "quick", compute_T = NULL)
#'summary_felogit(output)
#'
felogit <- function(data, formul = NULL, Option  = "quick", compute_X = "all", compute_T = "all", cluster = NULL, alpha = 0.05, CIOption = "CI2", nbCores = 4, varY = "Y", varX = NULL) {


  # Initialisation ----------------------------------------------------------
  nbDecimalsForRounding <- 4

  if (is.null(data)) { # Check data is there
    cat("Data should be provided.")
    output <- vector("list")
  } else {
    # Give default values if user set to NULL.
    # This is done for backward compatibility reasons.
    # Default algorithm is outer bounds. Option "quick" or "outer" both work (case insensitive)
    if (is.null(Option) || tolower(Option) %in% c("quick", "outer")) {
      Option <- "quick"
    }
    if (is.null(CIOption)) { # Default confidence interval is CI2
      CIOption <- "CI2"
    }
    if (is.null(cluster) && (any(class(data) != "environment") || !env_has(data, "clusterIndexes") || is.null(data$clusterIndexes))) { # Default is one observation per cluster
      if (is.null(formul)) {
        cluster <- 1:(dim(data$X)[1])
      } else {
        cluster <- 1:length(unique(data[[1]]))
      }
    } else if (is.null(cluster)) {
      cluster <- data$clusterIndexes
    }
    if (is.null(compute_T)) { # Default is to compute effects at all periods
      compute_T <- 1
    }
    if (!is.null(compute_X) && (length(compute_X) == 1)) {
      if (tolower(compute_X) == "all") {
        # Default is to compute effect for all variables
        compute_X <- NULL
      }
    }

  # Preprocess data ---------------------------------------------------------
    # Use formula, if specified, to turn data into wide format
    if (!is.null(formul)) {

      ### First column of data must be individual identifier
      varId <- colnames(data)[1]

      ### Second column of data must be periods
      varTime <- colnames(data)[2]

      ### Parse formula
      l <- attr(terms(formul), "variables")
      varY <- as.character(l[[2]])
      varX <- as.character(unlist(l[3:length(l)]))
      timeLabels <- sort(unique(data[[varTime]]))

      ### Turn into a wide array. We force ordering on individuals
      ### and periods to ensure we don't mix up data.
      nbIndls <- length(unique(data[[varId]]))
      nbPeriods <- length(unique(data[[varTime]]))
      dimX <- length(varX)
      X <- array(NA, c(nbIndls, nbPeriods, dimX))
      for (i in 1:dimX) {
        curVar <- varX[i]
        X[,, i] <- data[, c(varId, varTime, curVar)] |>
          pivot_wider(id_cols = {{varId}}, names_from = {{varTime}}, values_from = {{curVar}}) |>
          arrange({{varId}}) |>
          relocate(as.character(timeLabels)) |>
          select(!any_of(varId)) |>
          as.matrix()
      }
      Y <- data[, c(varId, varTime, varY)] |>
        pivot_wider(id_cols = {{varId}}, names_from = {{varTime}}, values_from = {{varY}}) |>
        arrange({{varY}}) |>
        relocate(as.character(timeLabels)) |>
        select(!any_of(varId)) |>
        as.matrix()
      colnames(Y) <- NULL
    } else { # If already in wide format, save variable names as specified in
      # the input, if specified

      Y <- data[[varY]]
      X <- data$X

      timeLabels <- colnames(data[[varY]])
      if (is.null(timeLabels)) {
        timeLabels <- 1:dim(data[[varY]])[2]
      }

      nbIndls <- dim(X)[1]
      if (length(dim(X)) == 2) {
        dimX <- 1
        X <- array(X, c(dim(X), dimX))
      } else {
        dimX <- dim(X)[3]
      }

      if (is.null(varX)) {
        if (is.null(dimnames(data$X)[[3]])) {
          varX <- paste0("X_", 1:(dim(data$X)[3]))
        } else {
          varX <- dimnames(data$X)[[3]]
        }
      }
    }

    # If one variable is not observed, all variables for the corresponding
    # individual-period pair are replaced by NAs.

    ### Find the degenerate individual-period pairs
    missingVariableAtIndlPeriodPair <- is.na(Y)
    for (i in 1:dimX) {
      missingVariableAtIndlPeriodPair <- missingVariableAtIndlPeriodPair | is.na(X[,, i])
    }

    ### Enforce NAs
    Y[missingVariableAtIndlPeriodPair] <- NA
    X[rep(c(missingVariableAtIndlPeriodPair), dimX)] <- NA


    # Discard constant variables
    varDiscard <- c()
    varKeep <- c()
    sdNA <- function(x) {return(sd(x, na.rm = TRUE))}

    for (i in 1:dimX) {
      ### Constant variables are those that don't vary across periods for over 99% of individuals
      if (quantile(apply(X[,, i], 1, sdNA), 0.99, na.rm = TRUE) == 0) {
        varDiscard <- c(varDiscard, i)
      } else {
        varKeep <- c(varKeep, i)
      }
    }

    ### Reshape X accordingly
    dimX <- length(varKeep)
    Tmax <- dim(X)[2]
    X <- array(X[,, varKeep], c(nbIndls, Tmax, dimX))

    ### Update labels and keep only the variables in compute_X that have not been discarded
    ### If compute_X is NULL, take all remaining variables
    discardedVarX <- varX[varDiscard]
    varX <- varX[varKeep]
    if (is.null(compute_X)) {
      compute_X <- varX
    } else if (is(compute_X,"numeric")) {
      compute_X <- varX[compute_X]
    } else {
      compute_X <- compute_X[compute_X %in% varX]
    }




    # Discard individuals with less than 2 observed periods
    indlsNotEnoughPeriods <- rowSums(!is.na(Y)) <= 1
    nbDiscardedIndls <- sum(indlsNotEnoughPeriods)
    Y <- matrix(Y[!indlsNotEnoughPeriods,], ncol = Tmax)
    X <- array(X[!indlsNotEnoughPeriods,,], dim(X) - c(nbDiscardedIndls, 0, 0))
    cluster <- cluster[!indlsNotEnoughPeriods]



    # Infer that variables where only two values are observed are binary
    continuousXvar <- rep(TRUE, dimX)
    for (i in 1:dimX) {
      observedValues <- unique(c(X[,, i]))
      observedValues <- observedValues[!is.na(observedValues)]
      continuousXvar[i] <- length(observedValues) > 2
    }

    ### If at least one of the variables for which we estimate the effect is
    ### binary, force the outer bound method
    if (any(varX[!continuousXvar] %in% compute_X)) {
      if (Option != "quick") {
        warning("One of the variables for which the effect is requested is binary. We reverted back to the outer bounds estimation method.")
        Option <- "quick"
      }
    }


    # Run estimation and format output ----------------------------------------
    estimators <- env()
    compute_X_bin <- which(varX %in% compute_X[compute_X %in% varX[!continuousXvar]])
    compute_X_cont <- which(varX %in% compute_X[compute_X %in% varX[continuousXvar]])
    data <- env("X" = X, "Y" = Y, "clusterIndexes" = cluster)
    if (length(compute_X_cont) != 0) {
      out_cont <- compute_AME(data, estimators, compute_X_cont, compute_T, Option, CIOption, alpha, nbCores)
    }
    if (length(compute_X_bin) != 0) {
      out_bin <- compute_ATE(data, estimators, compute_X_bin, compute_T, Option, CIOption, alpha, nbCores)
    }

    # Format CMLE results
    mat_results_CMLE <- matrix("", dimX, 4)
    colnames(mat_results_CMLE) <- c("Variable", "Point Est.", "se(Point Est.)", "pvalue")

    mat_results_CMLE[, 1] <- varX
    mat_results_CMLE[, 2] <- round(estimators$beta_hat$beta_hat, nbDecimalsForRounding)
    mat_results_CMLE[, 3] <- round(sqrt(diag(estimators$beta_hat$var_b) / (nbIndls - nbDiscardedIndls)), nbDecimalsForRounding)
    mat_results_CMLE[, 4] <- round(2 * (1 - pnorm(sqrt(nbIndls- nbDiscardedIndls) * abs(estimators$beta_hat$beta_hat) / sqrt(diag(estimators$beta_hat$var_b)))), nbDecimalsForRounding)

    # Format AME / ATE results
    # Look at the vignette to visualise the format we are building

    ### Give the right side to mat_results and write down each period
    ### on the first row of results relating to that period
    if (length(compute_T) == 0) {
      mat_results <- matrix("", length(compute_X), 5)
      mat_results[1, 1] <- "Tinf"
    } else if (compute_T[1] == "all") {
      mat_results <- matrix("", length(compute_X) * (Tmax + 1), 5)
      for (i in 1:Tmax) {
        mat_results[length(compute_X) * (i - 1) + 1, 1] <- paste0("T=", timeLabels[i])
      }
      mat_results[length(compute_X) * Tmax + 1, 1] <- "Average"
    } else {
      mat_results <- matrix("", length(compute_X) * length(compute_T), 5)
      for (i in 1:length(compute_T)) {
        mat_results[length(compute_X) * (i - 1) + 1, 1] <- paste0("T=", timeLabels[compute_T[i]])
      }
    }
    nbPeriods <- dim(mat_results)[1] / length(compute_X)

    ### Name the columns
    if (Option == "quick") {
      typeOfBounds <- "Estimated (outer) bounds"
    } else {
      typeOfBounds <- "Estimated (sharp) bounds"
    }
    colnames(mat_results) <- c("Period", "Variable", "AME/ATE", typeOfBounds, paste0(floor((1-alpha)*100), "% CI"))

    ### Name variables
    mat_results[, 2] <- rep(compute_X, nbPeriods)

    ### Decide whether AME or ATE
    mat_results[, 3] <- "AME"
    mat_results[rep(compute_X %in% varX[!continuousXvar], nbPeriods), 3] <- "ATE"

    ### Put in the estimated bounds and CI
    compute_X_bin <- which(compute_X %in% varX[compute_X_bin]) # Change the compute_X_...
    compute_X_cont <- which(compute_X %in% varX[compute_X_cont]) # to the rank inside compute_X

    for (i in 1:nbPeriods) {
      if (length(compute_X_bin) != 0) {
        mat_results[(i - 1) * length(compute_X) + compute_X_bin, 4] <- paste0("[", round(out_bin[[i]]$estimatedATEbounds[, 1], nbDecimalsForRounding), ", ", round(out_bin[[i]]$estimatedATEbounds[, 2], nbDecimalsForRounding), "]")
        mat_results[(i - 1) * length(compute_X) + compute_X_bin, 5] <- paste0("[", round(out_bin[[i]]$CI[, 1], nbDecimalsForRounding), ", ", round(out_bin[[i]]$CI[, 2], nbDecimalsForRounding), "]")
      }
      if (length(compute_X_cont) != 0) {
      mat_results[(i - 1) * length(compute_X) + compute_X_cont, 4] <- paste0("[", round(out_cont[[i]]$estimatedAMEbounds[, 1], nbDecimalsForRounding), ", ", round(out_cont[[i]]$estimatedAMEbounds[, 2], nbDecimalsForRounding), "]")
      mat_results[(i - 1) * length(compute_X) + compute_X_cont, 5] <- paste0("[", round(out_cont[[i]]$CI[, 1], nbDecimalsForRounding), ", ", round(out_cont[[i]]$CI[, 2], nbDecimalsForRounding), "]")
      }
    }

    # Write down the formula used as a string
    formul <- paste(varY, "~", do.call(paste, as.list(c(varX, "sep" = " + "))))

    if (!is.null(discardedVarX) && (length(discardedVarX) == 0)) {
      discardedVarX <- NULL
    }

    output <- vector("list")
    output[["summary"]] <- mat_results
    output[["n"]] <- nbIndls - nbDiscardedIndls
    output[["ndiscard"]] <- nbDiscardedIndls
    output[["Tmax"]] <- Tmax
    output[["vardiscard"]] <- discardedVarX
    output[["formul"]] <- formul
    output[["alpha"]] <- alpha
    output[["Option"]] <- Option
    output[["summary_CMLE"]] <- mat_results_CMLE
    output[["compute_T"]] <- compute_T
  }

  return(output)
}
