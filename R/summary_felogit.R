#' Produce the final summary table for the output of the felogit function
#'
#' @param output the output of the felogit function
#' @param format  can take value "latex" to print the latex table
#'
#' @return a kableExtra or xtable table plotted respectively in the R viewer or terminal
#' @export
#'
# @examples
summary_felogit <- function(output ,format = NULL){

  formula_p = paste0("Formula : ", as.character(output$formul)[[2]], " ~ ", as.character(output$formul)[[3]], ". \n ")
  discard_obs =  paste0("Nb of discarded individuals: ", output$ndiscard , ". \n ")
  obs =  paste0("Nb of observed individuals: ", output$n , ". \n ")
  maxT =  paste0("Maximal number of observed periods: ", output$Tmax , ". \n ")
  if(output$Option =="quick"){
    Option_c =  paste0("The method used to compute AME/ATE is the \"", output$Option, "\" method (i.e. the second method in DDL). \n ")
  }else{
    Option_c =  paste0("The method used to compute AME/ATE is the \"", output$Option, "\" method (i.e. the first method in DDL). \n ")
  }
  alpha_rep0 = paste0("The column period indicates at which period the AME/ATE are computed.  \n ")

   if(output$compute_T=="all"){
    alpha_rep0 = paste0(alpha_rep0, paste0("Average corresponds to the average of AME/ATE over the different periods. \n "))
  }

  alpha_rep =  paste0("CI: Confidence intervals are computed at a ", output$alpha*100  , "% level. \n ")
  general_note = paste0(formula_p, discard_obs,  obs,  maxT,  Option_c,  alpha_rep0 , alpha_rep)

  if(!is.null(output$vardiscard)){
    if(length(output$vardiscard)==1){
      Footnote_1 = paste0("The variable ", output$vardiscard, " has been discarded because it is constant or almost constant over time. \n")
    }else{
      Footnote_1 = paste0("The variables ", paste0(c(output$vardiscard), sep="",collapse = ", ") , " have been discarded because they are likely to be constant over time. \n")
    }
    general_note = paste0(general_note,Footnote_1 )
  }

  if(is.null(format)){

      cat("Estimates of coefficients in the fixed effect logit model (CMLE)")
      print(output$summary_CMLE %>%
                kable(caption = ,format='rst'))
      cat("Estimates of the Average Marginal Effects in the fixed effect logit model")
      print(output$summary %>%
        kable(caption = ,format='rst'))
       # kable_classic_2( html_font = "Cambria") %>%
      #  footnote(general = general_note))
      cat( general_note)


  }else if(format =="html"){

    general_note = paste0("Left table: Estimates of coefficients in the fixed effect logit model (CMLE) \n",
                          "Right table: Estimates of the Average Marginal Effects \n"
                          ,general_note)
    print(list(output$summary_CMLE,output$summary) %>%
            kable(caption = "Estimates of the Average Marginal Effects in the fixed effect logit model")  %>%
     kable_classic_2( html_font = "Cambria") %>%
      footnote(general = general_note))

    }else{

      general_note = paste0("Left table: Estimates of coefficients in the fixed effect logit model (CMLE) \n",
                            "Right table: Estimates of the Average Marginal Effects \n"
                            ,general_note)
      kable( list(output$summary_CMLE,output$summary), format = "latex", booktabs = T) %>%
        kable_styling() %>%
        footnote(general = general_note)
    }

  # cat("The table is printed in the viewer.")


}
