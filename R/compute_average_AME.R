#' Compute the average of the average marginal effects over all the periods between 1 and Tinf.
#'
#' @param output a list containing the outputs of the estimation of the bounds on the AME for the periods between 1 and Tinf.
#' @param Option chosen method to compute the AME: either sharp or quick. Default is quick.
#' @param Tinf vector of size n x1 containing the number of periods observed for each individual.
#' @param dimX the number of covariates
#' @param selectX the vector of selected covariates to compute the AME.
#' If NULL then bounds are computed for all covariates. NULL by default.
#' @param CIOption the option for the choice of the type of confidence intervals for the quick method, either CI2 or CI3. Default is CI2.
#' @param alpha the confidence level for the confidence intervals. Default is 5\%.
#' @param g_labels a matrix nx1 containing the individual labels refering to the type of attrition observed and described in the table G_types.
#' @param G_types a matrix describing the different possible types of attrition observed in the dataset
#' @param G_indic a matrix nxT containing the individual periods observed (0 if unobserved)
#'
#' @return  A list containing the values of : Option, n, Tmax,
#'
#'  - Delta_hat: either the estimator of the bounds on the average of Delta (sharp method) or the approximation of Delta (quick method) over the periods;
#'
#'  - length_CI: the length of the confidence intervals;
#'
#'  - et: the estimated standard error of the influence function(s), either of the two bounds on the average of Delta (sharp method) or of
#'  the approximation of Delta (quick method);
#'
#'  - bias_sup: the estimated upper bound on the bias (quick method);
#'
#'  - CI: the confidence intervals at the alpha level;
#'
#'  - b_hat: the estimated value of beta0 using the conditional maximum likelihood estimator;
#'
#'  - std_b: the estimated standard deviation of the inflence function of the  conditional maximum likelihood estimator of beta;
#'
#'  - influence: the matrix of size n x dimX containing the inflence function of the  conditional maximum likelihood estimator of beta;
#'
#' @export
#'
# @examples
compute_average_AME <- function(output,Option,Tinf,dimX,selectX,CIOption,alpha,g_labels , G_types, G_indic){

  ## find length of output.
  Tmax =  length(output)
  n = length(Tinf)

  ### number of variables to compute the AME/ATE
  if(is.null(selectX)){
    nb_var = dimX
    selectX = 1:dimX
  }else{
    nb_var = length(selectX)
  }


  # Option..
  if(Option=="quick"){

    # t_end=1

    Delta_hat = matrix(0, 1, dimX)
    bias_sup  = matrix(0, 1, dimX)
    infl = matrix(0,n,dimX)

    for(t_end in 1:Tmax){


      if(!is.null(G_types)){
        sel_g = G_types[,t_end]==1
        # discard observations if T not in periods.
        Tall = matrix(NA,dim( g_labels)[1],1)
        for(g in 1:length( sel_g)){
          if( sel_g[g]){
            Tall[g_labels==g,1]<-t_end
          }
        }
      }else{
        Tall = pmin(t_end,Tinf)
        # discard observations if T < t_end
        Tall[Tall < t_end] <- NA
      }



    out_cur = output[[t_end]]
    for_average = out_cur$for_average
    # Compute average of all estimators at period t
    Delta_hat =  Delta_hat + for_average[[1]]
    # Compute average of associated biases
    bias_sup   = bias_sup + for_average[[2]] # + out_cur$bias_sup
    # Compute average of influence functions
    # infl[Tinf>=t_end ,] = infl[Tinf>=t_end ,]  + out_cur$influence
    infl[!is.na(Tall),] = infl[!is.na(Tall) ,]  + out_cur$influence
    # dim(out_cur$influence)

    }
    # Delta_hat = Delta_hat/Tmax
    # bias_sup   =  bias_sup/Tmax
    infl = infl/(matrix(Tinf,n,1)%*%rep(1,dimX))
    ### Apply the IC construction method

    et = matrix(0,1,dimX)
    if(dimX==1){
      et =et + var(infl)
    }else{
      et =et + diag(var(infl))
    }

    et =sqrt(et)/sqrt(dim(infl)[1])
    et = et[selectX]
    # et = std(infl)/sqrt(dim(infl)[1])
    eps_n = (2*log(log(n)))^(1/2)/sqrt(dim(infl)[1])

    if( CIOption == "CI3"){
      # length_CI[s] = 2 * et[s] * sqrt(ncx2inv(0.95,1,(bias_sup[s]/et[s])^2));
      length_CI = 2 * et * sqrt(qchisq(0.95,1, ncp = ((bias_sup + eps_n)/et)^2, lower.tail = TRUE, log.p = FALSE))
    }else{
      length_CI = 2 * et * sqrt(qchisq(0.95,1, ncp = (bias_sup/et)^2, lower.tail = TRUE, log.p = FALSE))
    }
    ### compute CI

    CI =  matrix(Delta_hat, nb_var,1)%*%rep(1,2)
    CI[,1] = CI[,1] -  matrix(length_CI/2,  nb_var,1)
    CI[,2] = CI[,2] +  matrix(length_CI/2,  nb_var,1)

  }else{

    # t_end = 1

    Delta_hat = matrix(0, length(selectX),2)

    infl = vector("list")
    for(k in 1:dimX){
      infl[[k]]=matrix(0,n,2)
    }


    for(t_end in 1:Tmax){
      # attributes(out_cur)
      out_cur = output[[t_end]]
      for_average = out_cur$for_average
      # Compute average of all estimators at period t
      Delta_hat =  Delta_hat + for_average[[1]] #+ out_cur$Delta_hat
      # Compute average of influence functions
      for(k in 1:dimX){
        infl[[k]] = infl[[k]]  + out_cur$influence[[k]]
      }
      # dim( out_cur$influence[[k]])

    }

    for(k in 1:dimX){
      infl[[k]] = infl[[k]]/(matrix(Tinf,n,1)%*%rep(1,dim(infl[[k]])[2]))
    }
    # Delta_hat = Delta_hat/Tmax
    # infl = infl/(matrix(Tinf,n,1)%*%rep(1,dimX))
    ### Apply the IC construction method


    # Delta_hat[,1] - length_CI/2
    length_CI= matrix(NA,1, nb_var)
    et=matrix(NA,dimX,2)
    for(k in selectX){
      et[k,] =apply(infl[[k]],2,std)/sqrt(n);
    }
    et=matrix(et[selectX,], nb_var,2)
    bias_sup = matrix(0,2,1)

    time = 0;

    CI = matrix(NA, nb_var,2)

    ## compute CI at level alpha.
    # k=1
    for(k in 1:nb_var){
      quant = quantile_IM(alpha, Delta_hat[k,], et[k,]);

      # phi_alpha=1{sqrt(n)|beta hat_k|/std(phi_k)>q_{1-alpha/2}}
      phi_alpha = sqrt(n)*abs(out_cur$b_hat[k])/std(infl[[k]])> qnorm(1-alpha/2)
      if(phi_alpha){
        CI[k,] = cbind(Delta_hat[k,1] - quant* et[k,1], Delta_hat[k,2] + quant*et[k,2]);
      }else{
        CI[k,] = cbind(min(0,Delta_hat[k,1] - quant* et[k,1]), max(0,Delta_hat[k,2] + quant*et[k,2]));
      }

      # phi_alpha_m[k] = phi_alpha*1
      length_CI[k] = CI[k,2] - CI[k,1]

    }

    # std_b = apply(phi_b,2,std)/sqrt(dim(phi_b)[1])



  }


  out <- vector("list")
  out[[1]] <- Option
  out[[2]] <- n
  out[[3]] <- "average"
  out[[4]] <- NA
  out[[5]] <- Delta_hat
  out[[6]] <- length_CI
  out[[7]] <- et
  out[[8]] <- bias_sup
  out[[9]] <- CI
  out[[10]] <- NA
  out[[11]] <- out_cur$b_hat
  out[[12]] <- out_cur$std_b
  out[[13]] <- infl

  names(out) <- c("Option","n","Tmax","Time","Delta_hat","length_CI","et","bias_sup","CI","phi_alpha","b_hat","std_b", "influence")

 return(out)
}
