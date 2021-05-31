#' Function to perform the Monte-Carlo simulations displayed in DDL for the ATE using parallel computation.
#'
#' @param s indexes the  Monte-Carlo simulations
#' @param selectX the selected covariate with respect to which we compute the ATE
#' @param grid_T is the grid of values of maximal periods T observed in the simulated dataset
#' @param prop_T is the vector of proportions of invidividals with respective values of T in grid_T.
#' @param n is the number of individuals in the simulated dataset
#' @param dimX is the number of covariates
#' @param DGP is the DGP, with same numerotation as in DDL
#' @param beta0 is the vector of the true values of the parameter beta
#' @param Option the method choosen to perform estimation of the bounds
#' @param alpha_c_v the hypothetical vector of clustered values for alpha
#' @param nb_clust the number of hypothetical clusters
#' @param compute_T the vector of selected periods to compute the AME/ATE.
#' If NULL, then as described in Section 5.4 of DDL, ATE are computed at the mininum of the support of T.
#' If specified to ``all",  ATE are computed at all available periods but the average over the latter
#'  is also computed. Default is NULL.
#'
#' @return a list of vectors containing the different outcomes of the estimation (see compute_ATE).
#' @export
#'
# @examples
boot_ATE <- function(s, selectX=1, grid_T, prop_T, n,dimX,DGP, beta0, Option  = "quick", alpha_c_v =0 , nb_clust = 1,compute_T = NULL){
  # for( s in 1:nb_sim){
  #   if(s/100==round(s/100)){
  #     disp(s);
  #   }

  ## Generating the data
  # set.seed(102020)


  Tmax = max(grid_T)
  Xall =array(NA,c(n,Tmax,dimX))
  Yall = matrix(NA,n,Tmax)
  Tall = matrix(NA,n,1)
  # nt = floor(prop_T*n)
  # nindex = c(0,cumsum(nt))
  #
  nindex=0
  for(ci in 1:nb_clust){
    nindex = c(nindex,nindex[length(nindex)] + cumsum(floor(prop_T*n/nb_clust)))
  }
  nt = diff(nindex)

  Call = matrix(0,n,1)
  ## Generating the data
  # set.seed(102020)
  count = 1
  for(ci in 1:nb_clust){
    for(ti in 1:length(grid_T)){
      T = grid_T[ti]
      nt0 = nt[ti]

      #### stock the T_i
      Tall[(nindex[ti]+1):nindex[ti+1]] <- T

      X= (array(rand(nt0,Tmax*dimX),c(nt0,T,dimX))<=0.5)*1;

      # dim(  X)
      if( length(dim(X)) ==3){
        XT = matrix(X[,T,],nt0,dimX)
      }else{
        XT = matrix(X[,T],nt0,dimX)
      }

      Xbeta0 = zeros(nt0,T);

      # t= 1
      for( t in 1:T){
        Xbeta0[,t] = matrix(X[,t,],nt0,dimX)%*%beta0;
      }

      if(DGP==1){
        alpha = zeros(nt0,1)+ alpha_c_v[ci];
      }else if( DGP==2){
        alpha = XT[,1] + 2*(rand(nt0,1)<=.5)-1+ alpha_c_v[ci];
      }else if( DGP==3){
        alpha = XT[,1] + randn(nt0,1)+ alpha_c_v[ci];
      }else{
        prob_un = pmax(0,abs(XT[,1])-1/4);
        temp = rand(nt0,1);
        alpha = XT[,1] - (temp<=prob_un)*1 + (temp>1-prob_un)*1 + alpha_c_v[ci]
      }

      index=Xbeta0+ alpha%*%rep(1,dim(X)[2])
      prob = Lambda(index);
      Y = (rand(nt0,T)<=prob)*1;

      Xall[(nindex[count]+1):nindex[count+1],1:T, ] <-  X
      Yall[(nindex[count]+1):nindex[count+1], 1:T] <-  Y
      Call[(nindex[count]+1):nindex[count+1]] <- ci
      count =  count +1
    }
  }
  # Option  = "quick"
  # Option  = "continuous"
  # out = compute_AME(Y,X, start_point = start_point, Option  = Option , start_bounds= start_bounds)
  # case_DGP = 1
  out = compute_ATE(Yall,Xall, Call=Call, Option  = Option,selectX= selectX,compute_T= compute_T)

  # out[[4]] <- time
  # out[[5]] <- Delta_hat
  # out[[6]] <- length_CI
  # out[[7]] <- et
  # out[[8]] <- bias_sup
  return(out)

} # end of the for "s" simus
