#' Produces the final summary table for the output of the felogit function
#'
#' @param output the output of the felogit function
#' @param format can take value "latex" to print the latex table or "html" to
#' print as html
#'
#' @return a kableExtra or xtable table plotted respectively in the R viewer or terminal
#' @export
#'
# @examples
summary_felogit <- function(output, format = NULL) {

  formula_p <- paste0("Formula: ", output$formul, ".\n")
  discard_obs <- paste0("Number of discarded individuals: ", output$ndiscard , ".\n")
  obs <- paste0("Number of non-discarded individuals: ", output$n , ".\n")
  maxT <- paste0("Number of observed periods: ", output$Tmax , ".\n")
  if (output$Option == "quick") {
    Option_c <- paste0("The method used to compute the AME/ATE is the \"", output$Option, "\" method (i.e. the second method in DDL).\n")
  } else {
    Option_c <- paste0("The method used to compute AME/ATE is the \"", output$Option, "\" method (i.e. the first method in DDL).\n")
  }
  alpha_rep0 <- paste0("The column period indicates at which period the AME/ATE are computed.\n")

  if ((length(output$compute_T) == 1) && (output$compute_T == "all")) {
    alpha_rep0 <- paste0(alpha_rep0, paste0("Average corresponds to the average of AME/ATE over the different periods.\n"))
  }

  alpha_rep <- paste0("CI: Confidence intervals are computed at a ", round(output$alpha*100, 2), "% level.\n")
  general_note <- paste0(formula_p, discard_obs, obs, maxT, Option_c, alpha_rep0, alpha_rep)

  if (!is.null(output$vardiscard)) {
    if (length(output$vardiscard) == 1) {
      Footnote_1 <- paste0("The variable ", output$vardiscard, " has been discarded because it is constant or almost constant over time.\n")
    } else {
      Footnote_1 <- paste0("The variables ", paste0(c(output$vardiscard), sep = "", collapse = ", ") , " have been discarded because they are likely to be constant over time.\n")
    }
    general_note <- paste0(general_note, Footnote_1)
  }

  if (is.null(format)) {
  # Print in console
    cat("Estimates of coefficients in the fixed effect logit model (CMLE)")
    print(kable(output$summary_CMLE, caption = NULL, format = "rst"))
    cat("Estimates of the Average Marginal/Treatment Effects in the fixed effect logit model")
    print(kable(output$summary, caption = "", format = "rst"))
    cat(general_note)
  } else if (format == "html") {
  # Print in HTML
    general_note <- paste0(
      "Left table: Estimates of coefficients in the fixed effect logit model (CMLE)\n",
      "Right table: Estimates of the Average Marginal Effects\n",
      general_note
    )
    footnote(
      kable_classic_2(
        kable(
          list(output$summary_CMLE, output$summary),
          caption = "Estimates of the Average Marginal Effects in the fixed effect logit model",
          format = "html"
        ),
        html_font = "Cambria"
      ),
      general = general_note
    )
  } else {
  # Print in LaTeX
    general_note <- paste0(
      "Left table: Estimates of coefficients in the fixed effect logit model (CMLE)\n",
      "Right table: Estimates of the Average Marginal/Treatment Effects\n",
      general_note
    )
    footnote(
      kable_styling(
        kable(list(output$summary_CMLE, output$summary), format = "latex", booktabs = T)
      ),
      general = general_note
    )
  }
}
