#' Function which compute the ATE at different values of T according to the selected values of compute_T
#'
#' @param Yall matrix n x T containing the values of the dependent variable Yt
#' @param Xall array of dimensions n x T x dimX containing the values of the predictors at the different periods.
#' @param Call matrix n x 1 containing the identifiers of the clusters to which each individual belongs. Default is NULL
#' @param Option chosen method to compute the ATE: either sharp or quick. Default is quick.
#' @param selectX the vector of selected covariates to compute the AME.
#' If NULL then bounds are computed for all covariates. NULL by default.
#' @param compute_T the vector of selected periods to compute the ATE.
#' If NULL, then as described in Section 5.4 of DDL, ATE is computed at min supp (T).
#' If specified to ``all",  the ATE is computed at all available periods but the average over the latter is also computed. Default is NULL.
#' @param alpha  the confidence level for the confidence intervals. Default is 5\%.
#' @param CIOption the option for the choice of the type of confidence intervals for the quick method, either CI2 or CI3. Default is CI2.
#' @param g_labels a matrix nx1 containing the individual labels refering to the type of attrition observed and described in the table G_types.
#' @param G_types a matrix describing the different possible types of attrition observed in the dataset
#' @param G_indic a matrix nxT containing the individual periods observed (0 if unobserved)
#'
#' @return  a list of all the outputs of compute_ATE_t containing the different results of the estimation of the different values of T considered.
#' @export
#'
# @examples
compute_ATE <- function(Yall,Xall, Call= NULL,  Option  = "quick", selectX= 1,  compute_T = NULL,  alpha = 0.05, CIOption = "CI2",
                           g_labels = NULL , G_types = NULL, G_indic = NULL){

  ## Compute the Tinf = min( Supp(T)) for all individuals
  Tinf =  apply(Yall,1,isnot)
  Tmax = max(Tinf)

  ## Get the dimension of X
  if(length(dim(Xall))==3){
    dimX = dim(Xall)[3]
  }else{
    dimX=1
  }

  grid_T = NULL
  n_dist = NULL
  ### find the distibution of T in the population and sample size
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
    T = grid_T[t]
    M1 = repmat( matrix(seq(0,T)),1,T+1) -  repmat(seq(0,T),T+1,1)
    mat_C_k_T[[T]]  = choose(repmat(T-(0:T),T+1,1), M1)
    cheb_coeff[[T]]  = fliplr(t(coeff_Cheb01(T+1)));
  }


  options(warn=-1)
  b_lin = optim(par = rep(0,dimX) , lin_reg ,Y=Yall,X=Xall)$par
  start_point = b_lin/4
  options(warn=0)

  options(warn=-1)
  b_hat = optim(par = start_point, log_lik_FE_logit,Y=Yall,X=Xall)$par
  options(warn=0)


  # Compute the influence function of beta_hat. Useful for inference on
  # Delta, at the end.
  phi_b = infl_func_beta(b_hat,Yall, Xall, Call1);
  std_b = apply(phi_b,2,std)/sqrt(dim(phi_b)[1])



  # n <- dim(X)[1]
  # T <- dim(X)[2]
  # if(length(dim(X))==3){
  #   dimX = dim(X)[3]
  # }else{
  #   dimX=1
  # }
  #
  # Option = "slow"
  #  Option = "quick"
  output = vector("list")
  append_name <- function(x){return(paste0("T_",x))}
  ### compute only at last period
  if( length(compute_T)==0 ){
    Tall = Tinf
    output[[1]] <-  compute_ATE_t(Yall,Xall, prop_T,
                                           grid_T,n_dist,Tmax,Tall,Tinf,Call1,mat_C_k_T,
                                           cheb_coeff,b_hat,alpha, CIOption ,Option,dimX,selectX, phi_b, std_b)

    names(  output) <- c("Tinf")
  }else if( length(compute_T)>0 & compute_T[1]!="all"){
    ## compute for selected periods
    # t0=1
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

      output[[t0]] <-  compute_ATE_t(Yall,Xall, prop_T,
                                              grid_T,n_dist,Tmax,Tall,Tinf,Call1,mat_C_k_T,
                                              cheb_coeff,b_hat,alpha, CIOption ,Option,dimX,selectX, phi_b, std_b)

    }
    # output$T_1$Delta_hat
    # output$T_2$Delta_hat
    # output$T_3$Delta_hat

    names(  output) <- apply(matrix(compute_T,length(compute_T),1),1,append_name )
  }else{
    ## compute for all periods.

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
      output[[t_end]] <-  compute_ATE_t(Yall,Xall, prop_T,
                                                 grid_T,n_dist,Tmax,Tall,Tinf,Call1,mat_C_k_T,
                                                 cheb_coeff,b_hat,alpha, CIOption ,Option,dimX,selectX, phi_b, std_b)

    }

    names(  output) <- apply(matrix(1:Tmax,length(1:Tmax),1),1,append_name )

    ### add computation of average.

    output[[Tmax+1]]  <- compute_average_ATE(output,Option,Tinf,selectX,CIOption,alpha,g_labels, G_types , G_indic )

    names(  output) <- c(apply(matrix(1:Tmax,length(1:Tmax),1),1,append_name ),"average")

  }


  return( output)
}

















