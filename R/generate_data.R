#' Generate the simulated data in wide format according to the DGPs in DDL.
#'
#' @param grid_T a vector containing the different numbers of periods observed in the dataset.
#' @param prop_T a vector containin the proportions of the different numbers of periods observed in the dataset.
#' @param n the sample size.
#' @param dimX the number of covariates.
#' @param DGP the number indexing the chosen DGP, specified as in DDL.
#' @param beta0 a vector containing the true values of the parameter beta.
#' @param type_cont a vector of dummies coding for the different types of the covariates, either continous (1) or binary (0).
#' @param alpha_c_v a vector of length the number of specified clusters specifying the average value of alpha in these clusters.
#' @param nb_clust the number of hypothetical clusters.
#'
#' @return a list containing:
#'
#' - Xall: an array of size n x T x dimX containing the simulated values of the different covariates at all periods for all individuals;
#'
#' - Yall: an matrix of size n x T containing the simulated values of the dependent variable Y at all periods for all individuals;
#'
#' - Call: a vector of length n containing the index of the hypothetical clusters for all individuals. Default is a vector of ones.
#'
#' @export
#'
# @examples
generate_data <- function(grid_T, prop_T, n,dimX,DGP, beta0, type_cont, alpha_c_v =0 , nb_clust = 1){

  Tmax = max(grid_T)
  Xall =array(NA,c(n,Tmax,dimX))
  Yall = matrix(NA,n,Tmax)
  Tall = matrix(NA,n,1)
  nindex=0
  for(ci in 1:nb_clust){
    nindex = c(nindex,nindex[length(nindex)] + cumsum(floor(prop_T*n/nb_clust)))
  }
  nt = diff(nindex)
  # nindex = c(nindex,cumsum(nt))

  ## Generating the data
  # set.seed(102020)
  # ti=1

  Call = matrix(0,n,1)
  count = 1
  for(ci in 1:nb_clust){
    for(ti in 1:length(grid_T)){
      T = grid_T[ti]
      nt0 = nt[ti]

      #### stock the T_i
      Tall[(nindex[ti]+1):nindex[ti+1]] <- T


      nb_c = sum(type_cont)
      nb_b = dimX - nb_c
      if(nb_c >0){
        Xc= array(rand(nt0,Tmax*nb_c),c(nt0,T,nb_c))-.5
      }
      if(nb_b >0){
        Xb=(array(rand(nt0,Tmax*nb_b),c(nt0,T,nb_b))<=0.5)*1;
      }

      if(nb_c >0 &   nb_b >0){
        X= abind( Xc, Xb, along=3)
      }else if( nb_c >0 & nb_b ==0){
        X=  Xc
      }else{
        X=  Xb
      }





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
        alpha = zeros(nt0,1) + alpha_c_v[ci];
      }else if( DGP==2){
        alpha = XT[,1] + 2*(rand(nt0,1)<=.5)-1+ alpha_c_v[ci];
      }else if( DGP==3){
        alpha = XT[,1] + randn(nt0,1)+ alpha_c_v[ci];
      }else{
        prob_un = pmax(0,abs(XT[,1])-1/4);
        temp = rand(nt0,1);
        alpha = XT[,1] - (temp<=prob_un)*1 + (temp>1-prob_un)*1+ alpha_c_v[ci]
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

  output = vector("list")
  output[[1]] <- Xall
  output[[2]] <- Yall
  output[[3]] <- Call


  names(output) <- c("Xall","Yall","Call")
  return( output)
}
