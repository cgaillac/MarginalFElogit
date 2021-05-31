#' Function which compute the AME at different values of T according to the selected values of compute_T
#'
#'
#' @param Yall matrix n x T containing the values of the dependent variable Yt
#' @param Xall array of dimensions n x T x dimX containing the values of the predictors at the different periods.
#' @param Call matrix n x 1 containing the identifiers of the clusters to which each individual belongs. Default is NULL
#' @param Option chosen method to compute the AME: either sharp or quick. Default is quick.
#' @param selectX the vector of selected covariates to compute the AME.
#' If NULL then bounds are computed for all covariates. NULL by default.
#' @param compute_T the vector of selected periods to compute the AME.
#' If NULL, then as described in Section 5.4 of DDL, AME is computed at min supp (T).
#' If specified to ``all",  the AME is computed at all available periods but the average over the latter is also computed. Default is NULL.
#' @param alpha  the confidence level for the confidence intervals. Default is 5\%.
#' @param CIOption the option for the choice of the type of confidence intervals for the quick method, either CI2 or CI3. Default is CI2.
#' @param g_labels a matrix nx1 containing the individual labels referring to the type of attrition observed and described in the table G_types.
#' @param G_types a matrix describing the different possible types of attrition observed in the dataset
#' @param G_indic a matrix nxT containing the individual periods observed (0 if unobserved)
#' @param nbCores the number of cores used by the program to compute the AME for the ``sharp" method.
#' @param ratio the ratio R in DDL for the nonparametric estimator of the conditional moments of S
#'
#' @return a list of all the outputs of compute_AME_t containing the different results of the estimation of the different values of T considered.
#' @export
#'
# @examples
compute_AME<- function(Yall,Xall, Call= NULL, Option  = "quick", selectX = NULL, compute_T = NULL, alpha = 0.05, CIOption = "CI2",
                           g_labels = NULL , G_types = NULL, G_indic = NULL , nbCores = 4, ratio=10){
  # compute_T = NULL
  # selectX = NULL
  # alpha = 0.05
  # CIOption = "CI2"
  # nbCores = 4
  ## Compute the Tinf = min( Supp(T)) for all individuals
  Tinf =  apply(Yall,1,isnot)

  ## Max of the Tinf in dataset
  Tmax = max(Tinf)

  ## Get the dimension of X
  if(length(dim(Xall))==3){
    dimX = dim(Xall)[3]
  }else{
    dimX=1
  }

  ### find the distibution of Tinf in the population and sample size
  grid_T = NULL
  n_dist = NULL
  for(t in 1:Tmax){
    if(  sum(Tinf==t)>0){
      grid_T = c(  grid_T , t)
      n_dist = c( n_dist,sum(Tinf==t))
    }
  }
  prop_T = n_dist/sum(n_dist)

  ## number of clusters
  if(is.null(Call)){
    Call1 = matrix(1, dim(Yall)[1],1)
  }else{
    Call1 =  Call
  }

  ## compute combinatorial numbers at all possible dates in the dataset.
  mat_C_k_T= vector("list")
  cheb_coeff=vector("list")
  for(t in 1:length(grid_T)){
    T0 = grid_T[t]
    M1 = repmat( matrix(seq(0,T0)),1,T0+1) -  repmat(seq(0,T0),T0+1,1)
    mat_C_k_T[[T0]]  = choose(repmat(T0-(0:T0),T0+1,1), M1)
    cheb_coeff[[T0]]  = fliplr(t(coeff_Cheb01(T0+1)));
  }


  ## consider linear regression est. /4 as starting point
  options(warn=-1)
  b_lin = optim(par = rep(0,dimX) , lin_reg ,Y=Yall,X=Xall)$par
  start_point = b_lin/4
  options(warn=0)

  # if(dimX==1){
  #   b_hat = optimize( log_lik_FE_logit,start_bounds,Y=Yall,X=Xall )$minimum
  # }else{

  ### estimate loglikelihood.
  # ** insert catching errors ex: delete constant variables, ect..
  options(warn=-1)
  b_hat = optim(par = start_point, log_lik_FE_logit,Y=Yall,X=Xall)$par
  options(warn=0)
  # }

  # Compute the influence function of beta_hat. Useful for inference on
  # Delta, at the end.
  phi_b = infl_func_beta(b_hat,Yall, Xall, Call1);
  std_b = apply(phi_b,2,std)/sqrt(dim(phi_b)[1])

  ##
  output = vector("list")
  append_name <- function(x){return(paste0("T_",x))}
  ### compute only at last period
  if( length(compute_T)==0 ){
    Tall = Tinf
    output[[1]] <-  compute_AME_t(Yall,Xall, prop_T,
                                           grid_T,n_dist,Tmax,Tall,Tinf,Call1,mat_C_k_T,
                                           cheb_coeff,b_hat,alpha, CIOption ,Option,dimX,  selectX, phi_b, std_b , nbCores , ratio)

    names(  output) <- c("Tinf")
  }else if( length(compute_T)>0 & compute_T[1]!="all"){
    ## compute for selected periods
    # t0=1
    # compute_T= c(2,3)
    for(t0 in 1:length(compute_T)){

      t_end = compute_T[t0]
      ## Tall is the T at which the effect is computed, which is different according to the label g (if T is in the list of observed periods)
      # , g_labels  , G_types, G_indic)

      # find types containing t_end
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
      # cbind( Tall,dataX[,,1])
      grid_T0 = NULL
      n_dist0 = NULL
      for(t in 1:Tmax){
        if(  sum(Tinf[!is.na(Tall)]==t)>0){
          grid_T0 = c(  grid_T0 , t)
          n_dist0 = c( n_dist0,sum(Tinf[!is.na(Tall)]==t))
        }
      }
      prop_T0 = n_dist0/sum(n_dist0)

      # Tall = pmin(t_end,Tinf)
      # discard observations if T < t_end
      # Tall[Tall < t_end] <- NA

      # grid_T0 = grid_T =3
      # grid_T0[grid_T0>t_end] = t_end
      output[[t0]] <-  compute_AME_t(Yall,Xall, prop_T0,
                                              grid_T0,n_dist0,Tmax,Tall,Tinf,Call1,mat_C_k_T,
                                              cheb_coeff,b_hat,alpha, CIOption ,Option,dimX,  selectX , phi_b, std_b , nbCores, ratio)

    }
    # output$T_2$Delta_hat
    # output$T_3$Delta_hat

    names(  output) <- apply(matrix(compute_T,length(compute_T),1),1,append_name )
  }else{
    ## compute for all periods

    for(t_end in 1:Tmax){
      # find types containing t_end
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

      grid_T0 = NULL
      n_dist0 = NULL
      for(t in 1:Tmax){
        if(  sum(Tinf[!is.na(Tall)]==t)>0){
          grid_T0 = c(  grid_T0 , t)
          n_dist0 = c( n_dist0,sum(Tinf[!is.na(Tall)]==t))
        }
      }
      prop_T0 = n_dist0/sum(n_dist0)
      output[[t_end]] <-  compute_AME_t(Yall,Xall, prop_T0,
                                                 grid_T0,n_dist0,Tmax,Tall,Tinf,Call1,mat_C_k_T,
                                                 cheb_coeff,b_hat,alpha, CIOption ,Option,dimX,  selectX, phi_b, std_b  , nbCores, ratio)

    }

    ### add computation of average.

    output[[Tmax+1]]  <- compute_average_AME(output,Option,Tinf,dimX,selectX,CIOption,alpha,g_labels, G_types , G_indic )

    names(  output) <- c(apply(matrix(1:Tmax,length(1:Tmax),1),1,append_name ),"average")
    ## add name average


  }


  return( output)
}
