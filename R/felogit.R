#' This function implements the estimators of the bounds on the AME/ATE proposed in DDL.
#'
#' The function summary_felogit applied to the output of this felogit function prints the table containing the estimation results.
#'
#' @param data a data frame containing the panel data in long, i.e., one line for one individual-time observation.
#' The first column must be the individual identifier and the second column the temporal one. It can contains NA, which are treated as explained in the vignette.
#' @param formul an object of class "formula": a symbolic description of the model to be fitted.
#' Models for felogit are specified symbolically. A typical model has the form "response ~ terms`" where response is the (numeric) response vector
#' and terms is a series of terms which specifies a linear predictor for response. Note that predictors which are constant over time are discarded.
#' @param Option the method chosen to compute the AME/ATE, either "sharp" or "quick". Default is ``quick".
#' @param compute_X the vector of selected covariates to compute the AME/ATE. If NULL then bounds are computed for all covariates. NULL by default.
#' @param compute_T the vector of selected periods to compute the AME/ATE.
#' If NULL, then as described in Section 5.4 of DDL, AME/ATE are computed at the minimum of  supp(T). If specified to ``all",
#' AME/ATE are computed at all available periods but the average over the latter is also computed. Default is NULL.
#' @param cluster the vector of identifiers for the clusters if any.  Default is NULL.
#' @param alpha the confidence level for the confidence intervals. Default is 5\%.
#' @param CIOption the option for the choice of the type of confidence intervals for the "quick" method, either "CI2" or "CI3". Default is ``CI2".
#' @param nbCores the number of cores used by the program to compute the AME for the "sharp" method.
#' @param ratio the ratio R in DDL for the nonparametric estimator of the conditional moments of S
#' To reduce the computational time, this function can use several cores, in which case the library snowfall should be loaded first. By default, nbCores is set to 4.
#'
#'
#' @return A list containing:
#'
#'  - summary: a dataframe containing the estimation results,
#'
#'  - n: the number of used individuals,
#'
#'  - ndiscard: the number of discarded individuals,
#'
#'  - Tmax: the maximal number of periods of observed,
#'
#'  - vardiscard: the label of the discarded variables,
#'
#'  - formul: the formula used,
#'
#'  - alpha: the level used for the confidence intervals,
#'
#'  - Option: the method which is used.
#'
#'  - summary_CMLE : a dataframe containing the estimation results of the CMLE
#'
#' @export
#'
#' @examples
#'library(pglm)
#' data('UnionWage', package = 'pglm')
#'
#'UnionWage$union = (UnionWage$union =="yes")*1
#'UnionWage$rural = (UnionWage$rural =="yes")*1
#'UnionWage$black = (UnionWage$com =="black")*1 # used as test for discarded variable because constant
#'UnionWage$NorthEast =(UnionWage$region =="NorthEast")*1
#'sub <- UnionWage[UnionWage$year <1986,]
############## Estimate ##########################################
### formula
#'formul = as.formula(union ~  exper + married +black )
#'output <- felogit(formul , data= sub, Option  ="quick", compute_T =  NULL )
#'summary_felogit(output)
#'


felogit <- function(formul =NULL , data   , Option  = "quick", compute_X = NULL , compute_T = NULL, cluster = NULL , alpha = 0.05,
                    CIOption = "CI2", nbCores=4, ratio=10){


  if (is.null( formul )) {
    formul = NULL
  }

  if (is.null( data )) {
    cat("Data should be provided.")
    output = vector("list")
  }else{

    if (is.null(Option)) {
      Option  = "quick"
    }

    if (is.null( compute_X)) {
      compute_X = NULL
    }

    if (is.null( compute_T)) {
      compute_T = NULL
    }

    if (is.null(  cluster)) {
      cluster = NULL
    }
    if (is.null(alpha)) {
      alpha = 0.05
    }
    if (is.null( CIOption)) {
      CIOption = "CI2"
    }

    if (is.null( ratio)) {
      ratio = 10
    }


    if(!is.null(formul)){

    ### first column of data must be individual identifier
    var_id = colnames(data)[1]
    ### second column of data must be periods
    var_time =  colnames(data)[2]

    ### parse formula
    l <- as.list(attr(terms(formul), "variables"))[-1]
    var_y = as.character(l[[1]])
    var_x = matrix(1,1,length(l)-1)
    for(i in 1:(length(l)-1)){
      var_x[i] = as.character(l[[i+1]])
    }
    var_x =  c( var_x)
    dimX = length(var_x)

    ### check if non constant variables
    vardiscard=NULL

    indic = NULL
    #i=1
    sdna = function(x){return(sd(x,na.rm=TRUE))}

    for (i in 1:length(var_x)){
      test0 <- quantile(tapply(data[,var_x[i]], data[,var_id], FUN=sdna), 0.99, na.rm=T)
      if(test0 ==0){
        indic = c(indic, i)
      }
    }
    #i =4
    if(!is.null( indic)){
      vardiscard=var_x[indic]
      var_x = var_x[-c(indic)]
      dimX= length(var_x)
    }


    ### cleaning data and tests.
    ### transformation of the data to wide
    ### balanced panel with na
    s0<- make.pbalanced(data[,c(var_id,var_time,var_y, var_x)], balance.type="fill" )
    s0$sumNA <- (rowSums(is.na(s0[,c(var_x ,var_y)]))>0)*1

    # stocknb max period
    nbmax =  length(unique(s0[,var_time]))
    ##
    # ss <- as.data.frame(s0 %>% group_by("idcode") %>% sum("sumNA"))
    ss <- tapply(s0$sumNA, s0[,var_id], FUN=sum)
    ss = as.data.frame(cbind(rownames(ss),ss))
    colnames(ss ) <- c(var_id,"count_id")
    sub_1 = merge(s0[,c(var_id,var_time)],ss,by.x=var_id , by.y =var_id , all.x=TRUE)
    ## nb of unobserved period  <= nbmax-2 => nb of observed periods >=2
    if(max(as.numeric(sub_1$count_id))< (nbmax-1)){
      ## no attrition
      sub_select <- s0
      ndiscard = 0
    }else{
      indic = as.numeric(sub_1$count_id)>=(nbmax-1)
      sub_select <- s0[!indic,]
      sub_discard <- s0[indic,]
      # report number of discarded individuals
      ndiscard = length(unique(sub_discard[,var_id]))
    }



    # record all time labels
    times = sort(unique(sub_select[,var_time]))
    Tmax = length(times)

    ### second stage : reshape in wide
    sub_wide <- reshape(sub_select[,c(var_id,var_time,var_y, var_x)], idvar = var_id, timevar = var_time, direction = "wide")

    nobs = dim(sub_wide)[1]
    sub_wide[79,]

    labels.x = matrix(NA,dimX,Tmax)
    labels.y = matrix(NA,1,Tmax)
    for(i in 1:Tmax){
      labels.y[i] = paste0(var_y,".",times[i])
      for(j in 1:dimX){
        labels.x[j,i] = paste0(var_x[j],".",times[i])
      }
    }

    dataX = array(NA,c(dim(sub_wide)[1],Tmax,dimX))
    for(j in 1:dimX){
      dataX[,,j] = as.matrix(sub_wide[,labels.x[j,]],dim(sub_wide)[1],Tmax)
    }

    # dim(dataX)

    dataY = as.matrix(sub_wide[,labels.y],dim(sub_wide)[1],Tmax)
    ### handle attrition if any
    G_all = (!is.na(dataY))*1
    for(j in 1:dimX){
      G_all =  G_all*(!is.na(dataX[,,j]))*1
    }

    # dataY[79,]
    ## 0 = NA
    indic = !duplicated(G_all)
    G_types = matrix(G_all[ indic,], sum( indic), dim(G_all)[2])
    rownames(G_types) <- 1:dim(G_types)[1]
    g_labels = matrix(NA,dim(G_all)[1],1)
    # g=1
    for(g in 1:dim(G_types)[1]){
      g_labels [rowSums(G_all==(matrix(1,dim(G_all)[1],1)%*%matrix(G_types[g,],1, dim(G_types)[2] ))) == dim(G_all)[2]] <- g
    }

    ### transform data.
    Xall = array(NA,dim(dataX))
    Yall  = array(NA,dim(dataY))

    # g=3
    for(g in 1:dim(G_types)[1]){
      tinf= sum( G_types[g,]==1)
      Xall[g_labels==g,1:tinf,]   = dataX[g_labels==g, G_types[g,]==1,]
      Yall [g_labels==g,1:tinf]   = dataY[g_labels==g, G_types[g,]==1]
    }

    if(is.null(cluster)){
      Call = NULL
    }else{
      Call = sub_wide[,c(cluster)]
    }

    G_indic = G_types*(matrix(1,dim(G_types),1,1)%*%(1:Tmax))

  }else{ ## no formula, implicit that its organized and in wide.

    var_x=NULL
    ndiscard=0
    vardiscard=NULL


    ### transformation of the data to wide
    Yall = data$Yall
    Xall = data$Xall
    Call = data$Call

    nobs  =dim( Xall)[1]
    formul = vector("list")
    formul[[2]] = "Y"
    formul[[3]] = ""
    for(l in 1:dim(Xall)[3]){
      if(l==1){
        formul[[3]] = paste0(formul[[3]], "X",l)
      }else{
        formul[[3]] = paste0(formul[[3]], "+ X",l)
      }
    }

    if(!is.null(data$Gall)){
      g_labels <- data$Gall
    }else{
      g_labels <- NULL
    }

    if(!is.null(data$G_types)){
      G_types  <- data$G_types
      G_indic = G_types*(matrix(1,dim(G_types),1,1)%*%(1:Tmax))
    }else{
      G_types <- NULL
      G_indic <- NULL
    }


  }

  # G_types*(matrix(1,dim(G_types),1,1)%*%(1:Tmax))
  ## get the types of the variables
  Tinf =  apply(Yall,1,isnot)


  Tmax = max(Tinf)
  dimX <- dim(Xall)[3]

  ## default is continuous
  type_cont= matrix(1,dimX,1)
  # i=1
  for (i in 1:dimX){
    if(length(table(Xall[,,i]))==2){
      type_cont[i] = 0
    }
  }

  ref_c =  (1:dimX)[type_cont==1]
  ref_b =  (1:dimX)[type_cont==0]

  if(!is.null(var_x)){
    var_x_c = var_x[ref_c]
    var_x_b = var_x[ref_b]

  }else{
    var_x_c = NULL
    var_x_b = NULL
    if(length( ref_c)>0){
      for(j in 1:length( ref_c)){
        var_x_c = c(var_x_c, paste0("Xc",j))
      }
    }

    if(length( ref_b)>0){
      for(j in 1:length( ref_b)){
        var_x_b = c(var_x_b, paste0("Xb",j))
      }
    }
  }

  #### force to use quick method if at least one binary variable.

  if(sum(type_cont==0)>0){
    Option = "quick"
    #### add warning message
  }


  if(sum(type_cont)>0){
    out_c = compute_AME(Yall,Xall, Call= Call,  Option , selectX = ref_c, compute_T  , alpha , CIOption,  g_labels  , G_types, G_indic, nbCores,ratio )
  }
  out_b = vector("list")
  if(sum(type_cont==0)>0){
    for(i in ref_b){
      out_b[[i]] = compute_ATE(Yall,Xall, Call= Call, Option , selectX =i, compute_T , alpha , CIOption,  g_labels  , G_types, G_indic )
    }
  }



  mat_results_CMLE = matrix("",dimX,4);

  # sum(is.na(out_c$T_3$influence))
  ### compute at T_inf
  if( length(compute_T)==0 ){
    mat_results = matrix("",dimX,5);
    mat_results[1,1] = "Tinf"
  }else if( length(compute_T)>0 & compute_T[1]!="all"){
    ### compute at specified dates
    mat_results = matrix("",dimX*length(compute_T),5);
    for(i in 1:length(compute_T)){
      mat_results[dimX*(i-1)+1,1] = paste0("T=",compute_T[i])
    }
  }else{
    ### compute at all dates plus average
    ta = 1:Tmax
    t0 = length(1:Tmax) +1
    mat_results = matrix("",dimX*t0  ,5);
    for(i in 1:t0){
      if(i!=t0){
        mat_results[dimX*(i-1)+1,1] = paste0("T=",ta[i])
      }else{
        mat_results[dimX*(i-1)+1,1] = "Average"
      }
    }
  }

  ###
  if(Option =="quick"){

    count_c=1
    count_b=1
    ## for all the variables
    #i=1
    ind=1
    rnd= 4
    for(i in 1:dimX){
      # if continuous
      if(type_cont[i]==1){
        # out_c




        if( length(compute_T)==0 ){
          # names(  output) <- c("Tinf")
          ## insert est beta.
          mat_results[ind,2] = var_x_c[count_c]

          #### insert CMLE estimates
          mat_results_CMLE[i,1] = var_x_c[count_c]
          mat_results_CMLE[i,2] = round(out_c$Tinf$b_hat[i], rnd)
          mat_results_CMLE[i,3] =  round(out_c$Tinf$std_b[i], rnd)
          pval = 2*(1- pnorm(abs(out_c$Tinf$b_hat[i])/out_c$Tinf$std_b[i]))
          mat_results_CMLE[i,4] = round(pval, rnd)


          mat_results[ind,3] = "AME"
          # ind=  ind +1
          ## insert est Delta
          # mat_results[ind,2] = "\\underline{T}"
          mat_results[ind,4] =  round(out_c$Tinf$Delta_hat[count_c],  rnd)
          # mat_results[ind,7] =  round(out_c$Tinf$bias_sup[count_c],  rnd)
          mat_results[ind,5] =  paste0("[",round(out_c$Tinf$CI[count_c,1], rnd),",",round(out_c$Tinf$CI[count_c,2], rnd),"]")
          # mat_results[ind,9] =  round(out_c$Tinf$length_CI[count_c],  rnd)

          ind=  ind +1
          count_c = count_c+1
        }else if( length(compute_T)>0 & compute_T[1]!="all"){
          ## compute for selected periods
          # names(  output) <- apply(matrix(compute_T,length(compute_T),1),1,append_name )

          ind0 = ind
          # j=1

          out_cur = out_c[[1]]
          mat_results_CMLE[i,1] = var_x_c[count_c]
          mat_results_CMLE[i,2] = round(out_cur$b_hat[i], rnd)
          mat_results_CMLE[i,3] =  round(out_cur$std_b[i], rnd)
          pval = 2*(1- pnorm(abs(out_cur$b_hat[i])/out_cur$std_b[i]))
          mat_results_CMLE[i,4] = round(pval, rnd)


          for(j in 1:length(compute_T)){
            out_cur = out_c[[1]]
            mat_results[ind+dimX*(j-1),2] = var_x_c[count_c]
            # mat_results[ind+dimX*(j-1),3] = round(out_cur$b_hat[i], rnd)
            # mat_results[ind+dimX*(j-1),4] =  round(out_cur$std_b[i], rnd)
            # pval = 2*(1- pnorm(abs(out_cur$b_hat[i])/out_cur$std_b[i]))
            # mat_results[ind+dimX*(j-1),5] = round(pval, rnd)
            mat_results[ind+dimX*(j-1),3] = "AME"
            out_cur = out_c[[j]]
            ## insert est Delta
            # mat_results[ind,2] =  compute_T[j]
            mat_results[ind+dimX*(j-1),4] =  round(out_cur$Delta_hat[count_c],  rnd)
            # mat_results[ind,7] =  round(out_cur$bias_sup[count_c],  rnd)
            mat_results[ind+dimX*(j-1),5] =  paste0("[",round(out_cur$CI[count_c,1], rnd),",",round(out_cur$CI[count_c,2], rnd),"]")
            # mat_results[ind,9] =  round(out_cur$length_CI[count_c],  rnd)
          }
          ind=  ind +1
          count_c= count_c+1
        }else{
          ## compute for all periods.
          # for(t_end in 1:Tmax){
          # names(  output) <- apply(matrix(1:Tmax,length(1:Tmax),1),1,append_name )

          out_cur = out_c[[1]]
          mat_results_CMLE[i,1] = var_x_c[count_c]
          mat_results_CMLE[i,2] = round(out_cur$b_hat[i], rnd)
          mat_results_CMLE[i,3] =  round(out_cur$std_b[i], rnd)
          pval = 2*(1- pnorm(abs(out_cur$b_hat[i])/out_cur$std_b[i]))
          mat_results_CMLE[i,4] = round(pval, rnd)

          ind0 = ind
          for(j in 1:length(out_c)){
            out_cur = out_c[[1]]
            mat_results[ind+dimX*(j-1),2] = var_x_c[count_c]
            # mat_results[ind+dimX*(j-1),3] = round(out_cur$b_hat[i], rnd)
            # mat_results[ind+dimX*(j-1),4] =  round(out_cur$std_b[i], rnd)
            # pval = 2*(1- pnorm(abs(out_cur$b_hat[i])/out_cur$std_b[i]))
            # mat_results[ind+dimX*(j-1),5] = round(pval, rnd)
            mat_results[ind+dimX*(j-1),3] = "AME"

            out_cur = out_c[[j]]
            ## insert est Delta
            # mat_results[ind,2] =  compute_T[j]
            mat_results[ind+dimX*(j-1),4] =  round(out_cur$Delta_hat[count_c],  rnd)
            # mat_results[ind,7] =  round(out_cur$bias_sup[count_c],  rnd)
            mat_results[ind+dimX*(j-1),5] =  paste0("[",round(out_cur$CI[count_c,1], rnd),",",round(out_cur$CI[count_c,2], rnd),"]")
            # mat_results[ind,9] =  round(out_cur$length_CI[count_c],  rnd)
          }
          ind=  ind +1
          count_c= count_c+1



        }

      }else{ # end type continous variable


        out_cur =  out_b[[i]]


        if( length(compute_T)==0 ){
          # names(  output) <- c("Tinf")
          ## insert est beta.
          mat_results[ind,2] = var_x_b[count_b]

          ### insert results CMLE
          mat_results_CMLE[i,1] = var_x_b[count_b]
          mat_results_CMLE[i,2] = round(out_cur$Tinf$b_hat[i], rnd)
          mat_results_CMLE[i,3] =  round(out_cur$Tinf$std_b[i], rnd)
          pval = 2*(1- pnorm(abs(out_cur$Tinf$b_hat[i])/out_cur$Tinf$std_b[i]))
          mat_results_CMLE[i,4] = round(pval, rnd)

          mat_results[ind,3] = "ATE"
          ## insert est Delta
          # mat_results[ind,2] = "\\underline{T}"
          mat_results[ind,4] =  round(out_cur$Tinf$Delta_hat,  rnd)
          # mat_results[ind,7] =  round(out_cur$Tinf$bias_sup,  rnd)
          mat_results[ind,5] =  paste0("[",round(out_cur$Tinf$CI[,1], rnd),",",round(out_cur$Tinf$CI[,2], rnd),"]")
          # mat_results[ind,9] =  round(out_cur$Tinf$length_CI,  rnd)
          ind=  ind +1
          count_b= count_b+1

        }else if( length(compute_T)>0 & compute_T[1]!="all"){
          ## compute for selected periods
          # names(  output) <- apply(matrix(compute_T,length(compute_T),1),1,append_name )

          out_cur0 = out_cur[[1]]
          ###
          mat_results_CMLE[i,1] = var_x_b[count_b]
          mat_results_CMLE[i,3] = round(out_cur0 $b_hat[i], rnd)
          mat_results_CMLE[i,3] =  round(out_cur0 $std_b[i], rnd)
          pval = 2*(1- pnorm(abs(out_cur0 $b_hat[i])/out_cur0 $std_b[i]))
          mat_results_CMLE[i,4] = round(pval, rnd)
          # ind=  ind +1


          ind0 = ind
          for(j in 1:length(compute_T)){
            mat_results[ind+dimX*(j-1),2] = var_x_b[count_b]
            mat_results[ind+dimX*(j-1),3] = "ATE"

            out_cur0 = out_cur[[j]]
            ## insert est Delta
            # mat_results[ind,1] =  compute_T[j]
            mat_results[ind+dimX*(j-1),4] =  round(out_cur0$Delta_hat,  rnd)
            # mat_results[ind,7] =  round(out_cur0$bias_sup,  rnd)
            mat_results[ind+dimX*(j-1),5] =  paste0("[",round(out_cur0$CI[,1], rnd),",",round(out_cur0$CI[,2], rnd),"]")
            # mat_results[ind,9] =  round(out_cur0$length_CI,  rnd)

            # }
          }
          ind=  ind +1
          count_b= count_b+1
        }else{
          ## compute for all periods.
          # for(t_end in 1:Tmax){
          # names(  output) <- apply(matrix(1:Tmax,length(1:Tmax),1),1,append_name )


          out_cur0 = out_cur[[1]]
          mat_results_CMLE[i,1] = var_x_b[count_b]
          mat_results_CMLE[i,2] = round(out_cur0 $b_hat[i], rnd)
          mat_results_CMLE[i,3] =  round(out_cur0 $std_b[i], rnd)
          pval = 2*(1- pnorm(abs(out_cur0 $b_hat[i])/out_cur0 $std_b[i]))
          mat_results_CMLE[i,4] = round(pval, rnd)

          ind0 = ind
          for(j in 1:length(out_cur)){
            mat_results[ind+dimX*(j-1),2] = var_x_b[count_b]
            # mat_results[ind+dimX*(j-1),3] = round(out_cur0 $b_hat[i], rnd)
            # mat_results[ind+dimX*(j-1),4] =  round(out_cur0 $std_b[i], rnd)
            # pval = 2*(1- pnorm(abs(out_cur0 $b_hat[i])/out_cur0 $std_b[i]))
            # mat_results[ind+dimX*(j-1),5] = round(pval, rnd)
            mat_results[ind+dimX*(j-1),3] = "ATE"
            # ind=  ind +1
            out_cur0 = out_cur[[j]]
            ## insert est Delta
            # mat_results[ind,1] =  compute_T[j]
            mat_results[ind+dimX*(j-1),4] =  round(out_cur0$Delta_hat,  rnd)
            # mat_results[ind,7] =  round(out_cur0$bias_sup,  rnd)
            mat_results[ind+dimX*(j-1),5] =  paste0("[",round(out_cur0$CI[,1], rnd),",",round(out_cur0$CI[,2], rnd),"]")
            # mat_results[ind,9] =  round(out_cur0$length_CI,  rnd)

            # }
          }
          ind=  ind +1
          count_b= count_b+1


        }




      }
    } # end dimX

    colnames(mat_results) <- c("Period","Variable","AME/ATE", "Estimate", paste0((1-alpha)*100,"% CI"))
    # colnames(mat_results) <- c("","Variable","\\widehat{\\beta}_k","\\widehat{\\sigma}", "","\\widehat{\\Delta}_k", "Bias(\\widehat{\\Delta}_k)", "CI_{\\alpha}", "CI length")

  }else{


    #
    # if( length(compute_T)==0 ){
    #   mat_results = matrix("",dimX*(1+1),8);
    # }else if( length(compute_T)>0 & compute_T[1]!="all"){
    #   mat_results = matrix("",dimX*(1+length(compute_T)),8);
    # }else{
    #   mat_results = matrix("",dimX*(1+lenght(1:Tmax)),8);
    # }
    #

    count_c=1
    count_b=1
    ## for all the variables
    #i=1
    ind=1
    rnd= 4
    for(i in 1:dimX){
      # if continuous
      if(type_cont[i]==1){
        # out_c
        if( length(compute_T)==0 ){
          # names(  output) <- c("Tinf")
          ## insert est beta.
          mat_results[ind,2] = var_x_c[count_c]

          mat_results_CMLE[i,1] = var_x_c[count_c]
          mat_results_CMLE[i,2] = round(out_c$Tinf$b_hat[i], rnd)
          mat_results_CMLE[i,3] =  round(out_c$Tinf$std_b[i], rnd)
          pval = 2*(1- pnorm(abs(out_c$Tinf$b_hat[i])/out_c$Tinf$std_b[i]))
          mat_results_CMLE[i,4] = round(pval, rnd)
          # ind=  ind +1

          mat_results[ind,3] = "AME"
          ## insert est Delta
          # mat_results[ind,2] = "\\underline{T}"
          # mat_results[ind,6] =  round(out_c$Tinf$Delta_hat,  rnd)
          mat_results[ind,4] =  paste0("[",round(out_c$Tinf$Delta_hat[count_c,1], rnd),",",round(out_c$Tinf$Delta_hat[count_c,2], rnd),"]")
          mat_results[ind,5] =  paste0("[",round(out_c$Tinf$CI[count_c,1], rnd),",",round(out_c$Tinf$CI[count_c,2], rnd),"]")
          # mat_results[ind,8] =  round(out_c$Tinf$length_CI[count_c],  rnd)
          ind=  ind +1
          count_c= count_c+1
        }else if( length(compute_T)>0 & compute_T[1]!="all"){
          ## compute for selected periods
          # names(  output) <- apply(matrix(compute_T,length(compute_T),1),1,append_name )

          out_cur = out_c[[1]]
          mat_results_CMLE[i,1] = var_x_c[count_c]
          mat_results_CMLE[i,2] = round(out_cur$b_hat[i], rnd)
          mat_results_CMLE[i,3] =  round(out_cur$std_b[i], rnd)
          pval = 2*(1- pnorm(abs(out_cur$b_hat[i])/out_cur$std_b[i]))
          mat_results_CMLE[i,4] = round(pval, rnd)

          for(j in 1:length(compute_T)){
            out_cur = out_c[[1]]
            mat_results[ind+dimX*(j-1),2] = var_x_c[count_c]
            # mat_results[ind+dimX*(j-1),3] = round(out_cur$b_hat[i], rnd)
            # mat_results[ind+dimX*(j-1),4] =  round(out_cur$std_b[i], rnd)
            # pval = 2*(1- pnorm(abs(out_cur$b_hat[i])/out_cur$std_b[i]))
            # mat_results[ind+dimX*(j-1),5] = round(pval, rnd)
            mat_results[ind+dimX*(j-1),3] = "AME"
            out_cur = out_c[[j]]
            ## insert est Delta
            # mat_results[ind,2] =  compute_T[j]
            # mat_results[ind,6] =  round(out_cur$Delta_hat,  rnd)
            mat_results[ind+dimX*(j-1),4] =  paste0("[",round(out_cur$Delta_hat[count_c,1], rnd),",",round(out_cur$Delta_hat[count_c,2], rnd),"]")
            mat_results[ind+dimX*(j-1),5] =  paste0("[",round(out_cur$CI[count_c,1], rnd),",",round(out_cur$CI[count_c,2], rnd),"]")
            # mat_results[ind,8] =  round(out_cur$length_CI[count_c],  rnd)
          }
          ind=  ind +1
          count_c= count_c+1
        }else{
          ## compute for all periods.
          # for(t_end in 1:Tmax){
          # names(  output) <- apply(matrix(1:Tmax,length(1:Tmax),1),1,append_name )

          out_cur = out_c[[1]]
          mat_results_CMLE[i,1] = var_x_c[count_c]
          mat_results_CMLE[i,2] = round(out_cur$b_hat[i], rnd)
          mat_results_CMLE[i,3] =  round(out_cur$std_b[i], rnd)
          pval = 2*(1- pnorm(abs(out_cur$b_hat[i])/out_cur$std_b[i]))
          mat_results_CMLE[i,4] = round(pval, rnd)

          for(j in 1:length(out_c)){
            out_cur = out_c[[1]]
            mat_results[ind+dimX*(j-1),2] = var_x_c[count_c]
            # mat_results[ind+dimX*(j-1),3] = round(out_cur$b_hat[i], rnd)
            # mat_results[ind+dimX*(j-1),4] =  round(out_cur$std_b[i], rnd)
            # pval = 2*(1- pnorm(abs(out_cur$b_hat[i])/out_cur$std_b[i]))
            # mat_results[ind+dimX*(j-1),5] = round(pval, rnd)
            mat_results[ind+dimX*(j-1),3] = "AME"
            out_cur = out_c[[j]]
            ## insert est Delta
            # mat_results[ind,2] =  compute_T[j]
            # mat_results[ind,6] =  round(out_cur$Delta_hat,  rnd)
            mat_results[ind+dimX*(j-1),4] =  paste0("[",round(out_cur$Delta_hat[count_c,1], rnd),",",round(out_cur$Delta_hat[count_c,2], rnd),"]")
            mat_results[ind+dimX*(j-1),5] =  paste0("[",round(out_cur$CI[count_c,1], rnd),",",round(out_cur$CI[count_c,2], rnd),"]")
            # mat_results[ind,8] =  round(out_cur$length_CI[count_c],  rnd)
          }
          ind=  ind +1
          count_c= count_c+1

        }

      }else{ # end type cont


        out_cur =  out_b[[i]]
        mat_results_CMLE[ind,1] =var_x_b[count_b]
        mat_results_CMLE[ind,2] = round(out_cur$Tinf$b_hat[i], rnd)
        mat_results_CMLE[ind,3] =  round(out_cur$Tinf$std_b[i], rnd)
        pval = 2*(1- pnorm(abs(out_cur$Tinf$b_hat[i])/out_cur$Tinf$std_b[i]))
        mat_results_CMLE[ind,4] = round(pval, rnd)

        if( length(compute_T)==0 ){
          # names(  output) <- c("Tinf")
          ## insert est beta.
          mat_results[ind,2] =var_x_b[count_b]
          # mat_results[ind,3] = round(out_cur$Tinf$b_hat[i], rnd)
          # mat_results[ind,4] =  round(out_cur$Tinf$std_b[i], rnd)
          # pval = 2*(1- pnorm(abs(out_cur$Tinf$b_hat[i])/out_cur$Tinf$std_b[i]))
          # mat_results[ind,5] = round(pval, rnd)
          mat_results[ind,3] = "ATE"
          # ind=  ind +1

          ## insert est Delta
          # mat_results[ind,2] = "\\underline{T}"
          # mat_results[ind,6] =  round(out_cur$Tinf$Delta_hat,  rnd)
          mat_results[ind,4] =  paste0("[",round(out_cur$Tinf$Delta_hat[,1], rnd),",",round(out_cur$Tinf$Delta_hat[,2], rnd),"]")
          mat_results[ind,5] =  paste0("[",round(out_cur$Tinf$CI[,1], rnd),",",round(out_cur$Tinf$CI[,2], rnd),"]")
          # mat_results[ind,8] =  round(out_cur$Tinf$length_CI,  rnd)
          ind=  ind +1
          count_b= count_b+1
        }else if( length(compute_T)>0 & compute_T[1]!="all"){
          ## compute for selected periods
          # names(  output) <- apply(matrix(compute_T,length(compute_T),1),1,append_name )


          out_cur0 = out_cur[[1]]
          mat_results_CMLE[i,1] =var_x_b[count_b]
          mat_results_CMLE[i,2] = round(out_cur0$b_hat[i], rnd)
          mat_results_CMLE[i,3] =  round(out_cur0$std_b[i], rnd)
          pval = 2*(1- pnorm(abs(out_cur0$b_hat[i])/out_cur0$std_b[i]))
          mat_results_CMLE[i,4] = round(pval, rnd)

          for(j in 1:length(compute_T)){
            mat_results[ind+dimX*(j-1),2] =var_x_b[count_b]
            # mat_results[ind+dimX*(j-1),3] = round(out_cur0$b_hat[i], rnd)
            # mat_results[ind+dimX*(j-1),4] =  round(out_cur0$std_b[i], rnd)
            # pval = 2*(1- pnorm(abs(out_cur0$b_hat[i])/out_cur0$std_b[i]))
            # mat_results[ind+dimX*(j-1),5] = round(pval, rnd)
            mat_results[ind+dimX*(j-1),3] = "ATE"
            # ind=  ind +1

            # for(j in 1:length(compute_T)){
            out_cur0 = out_cur[[j]]
            ## insert est Delta
            # mat_results[ind,1] =  compute_T[j]
            # mat_results[ind,6] =  round(out_cur0$Delta_hat,  rnd)
            mat_results[ind+dimX*(j-1),4] =  paste0("[",round(out_cur0$Delta_hat[,1], rnd),",",round(out_cur0$Delta_hat[,2], rnd),"]")
            mat_results[ind+dimX*(j-1),5] =  paste0("[",round(out_cur0$CI[,1], rnd),",",round(out_cur0$CI[,2], rnd),"]")
            # mat_results[ind,8] =  round(out_cur0$length_CI,  rnd)
          }
          ind=  ind +1
          count_b= count_b+1

        }else{
          ## compute for all periods.
          # for(t_end in 1:Tmax){
          # names(  output) <- apply(matrix(1:Tmax,length(1:Tmax),1),1,append_name )


          out_cur0 = out_cur[[1]]
          mat_results_CMLE[i,1] =var_x_b[count_b]
          mat_results_CMLE[i,2] = round(out_cur0$b_hat[i], rnd)
          mat_results_CMLE[i,3] =  round(out_cur0$std_b[i], rnd)
          pval = 2*(1- pnorm(abs(out_cur0$b_hat[i])/out_cur0$std_b[i]))
          mat_results_CMLE[i,4] = round(pval, rnd)

          for(j in 1:length(out_cur)){
            mat_results[ind+dimX*(j-1),2] =var_x_b[count_b]
            # mat_results[ind+dimX*(j-1),3] = round(out_cur0$b_hat[i], rnd)
            # mat_results[ind+dimX*(j-1),4] =  round(out_cur0$std_b[i], rnd)
            # pval = 2*(1- pnorm(abs(out_cur0$b_hat[i])/out_cur0$std_b[i]))
            # mat_results[ind+dimX*(j-1),5] = round(pval, rnd)
            mat_results[ind+dimX*(j-1),3] = "ATE"
            # ind=  ind +1

            # for(j in 1:length(compute_T)){
            out_cur0 = out_cur[[j]]
            ## insert est Delta
            # mat_results[ind,1] =  compute_T[j]
            # mat_results[ind,6] =  round(out_cur0$Delta_hat,  rnd)
            mat_results[ind+dimX*(j-1),4] =  paste0("[",round(out_cur0$Delta_hat[,1], rnd),",",round(out_cur0$Delta_hat[,2], rnd),"]")
            mat_results[ind+dimX*(j-1),5] =  paste0("[",round(out_cur0$CI[,1], rnd),",",round(out_cur0$CI[,2], rnd),"]")
            # mat_results[ind,8] =  round(out_cur0$length_CI,  rnd)
          }
          ind=  ind +1
          count_b= count_b+1


        }

      }
    } # end dimX





    colnames(mat_results) <- c("Period","Variable","AME/ATE", "Estimate", paste0((1-alpha)*100,"% CI"))
    # colnames(mat_results) <- c("Variable","T","\\widehat{\\beta}_k","\\widehat{\\sigma}", "","[\\underline{\\widehat{\\Delta}}_k, \\overline{\\widehat{\\Delta}}_k]", "CI_{\\alpha}", "CI length")



  }# end option quick



  colnames(mat_results_CMLE) <- c("Variable","Point Est.","se(Point Est.)", "pvalue")
  ### format the output for the print function.

  output = vector("list")
  output[[1]] <- mat_results
  output[[2]] <- nobs
  output[[3]] <- ndiscard
  output[[4]] <- Tmax
  output[[5]] <- vardiscard
  output[[6]] <- formul
  output[[7]] <- alpha
  output[[8]] <- Option
  output[[9]] <- mat_results_CMLE
  output[[10]] <- if(is.null(compute_T)){"NULL"}else{compute_T}
  names(output) <- c("summary", "n", "ndiscard", "Tmax", "vardiscard", "formul", "alpha", "Option","summary_CMLE", "compute_T")
  }

  return(output)
}

















