#' Compute the influence function of beta hat: E(d^2ln L/dbeta^2)^{-1}  x (ln L/dbeta)(Yi|Xi,beta).
#'
#' @param beta the parameter at which the influence function is computed
#' @param Y a matrix of size nxT containing the values of the dependent variable Y
#' @param X an array of size nxTxdimX containing the values of the covariates X
#' @param Call a vector of size n containing the indexes of the hypothetical clusters for the individuals
#'
#' @return a vector of size nx1 containing the value of the influence function
#' @export
#'
# @examples
infl_func_beta <- function(beta,Y, X, Call=NULL){


  # [n, T, k]=size(X);
  ## change for dim(X) >1
  n <- dim(X)[1]
  Tmax <- dim(X)[2]
  # Tmax = T
  if(length(dim(X))==3){
    k = dim(X)[3]
  }else{
    k=1
  }
  S = apply(Y,1,sumNA);
  Tall =  apply(Y,1,isnot)
  grid_T = NULL
  ### find the distibution of T in the population and sample size
  for(t in 1:Tmax){
    if(  sum(Tall==t)>0){
      grid_T = c(  grid_T , t)
    }
  }


  V = zeros(n,Tmax);
  # if(k==1){
  #   V = exp(X[,,1]*beta);
  # }else{
  #   for( t in 1:T){
  #     V[,t] =  exp(X[,t,]%*%matrix(beta));
  #   }
  # }

  index = zeros(n,Tmax)
  ## test if X of dim > 2
  if(length(dim(X))==3){
    for(ti in grid_T){
      for (t in 1:ti){
        index[Tall==ti,t] = X[Tall==ti,t,]%*%matrix(beta);
      }
    }
  }else{
    index = X*beta;
  }
  V = exp(index);

  score = zeros(k,n);
  # hess = zeros(k);


  denom = matrix(0,n,1)
  for(ti in grid_T){
    denom[Tall==ti,1] = C_S_fun(S[Tall==ti],matrix(V[Tall==ti,1:ti], sum(Tall==ti) , ti));
  }

  J_denom = matrix(NA,n,Tmax);
  hess_denom = array(NA,c(n,Tmax,Tmax));

  XX = array(NA,c(n,Tmax,Tmax));


  if(is.null(Call)){

  for(ti in grid_T){
    if(ti==1){
      t= 1
      J_denom[Tall==ti,t] = V[Tall==ti,t]/denom[Tall==ti];
    }else{
    for( t in seq(1,ti)){
      V0 = matrix(V[Tall==ti& Call ==ci,1:ti],sum(Tall==ti& Call ==ci),ti)
      J_denom[Tall==ti,t] = V[Tall==ti,t]*C_S_fun(S[Tall==ti]-1, matrix(V0[,(1:ti)!=t],sum(Tall==ti),sum((1:ti)!=t)) )/denom[Tall==ti];
    }
    }
  }
  # % Computation of (ln L/dbeta)(Yi|Xi,beta).

  if(k==1){
    score[1,] = t(apply(X[,,1]*(Y - J_denom),1,sumNA));
  }else{
    for(j in 1:k){
      score[j,] = t(apply(X[,,j]*(Y - J_denom),1,sumNA));
    }
  }

  # for( i in seq(1,k)){
  # score[i,] = t(apply(X[,,i]*(Y - J_denom),1,sum));

  # }

  info_Fish = cov(t(score));

  }else{

    nb_clust =length(unique(Call))
    info_Fish = 0
    for(ci in 1: nb_clust){
      for(ti in grid_T){
        if(ti==1){
          J_denom[Tall==ti & Call ==ci,t] = V[Tall==ti& Call ==ci,t]/denom[Tall==ti& Call ==ci];
        }else{
          for( t in seq(1,ti)){
            V0 = matrix(V[Tall==ti& Call ==ci,1:ti],sum(Tall==ti& Call ==ci),ti)
            J_denom[Tall==ti & Call ==ci,t] = V[Tall==ti& Call ==ci,t]*C_S_fun(S[Tall==ti& Call ==ci]-1, matrix(V0[,(1:ti)!=t],sum(Tall==ti& Call ==ci),sum((1:ti)!=t)) )/denom[Tall==ti& Call ==ci];
          }
        }
      }



    # % Computation of (ln L/dbeta)(Yi|Xi,beta).

    if(k==1){
      score[1,Call ==ci] = t(apply(matrix(X[Call ==ci,,1],sum(Call ==ci) , dim(X)[2] )*(Y[Call ==ci,] - J_denom[Call ==ci,]),1,sumNA));
    }else{
      for(j in 1:k){
        score[j,Call ==ci] = t(apply(matrix(X[Call ==ci,,j],sum(Call ==ci) ,  dim(X)[2] )*(Y[Call ==ci,] - J_denom[Call ==ci,]),1,sumNA));
      }
    }

    info_Fish = info_Fish + cov(t(score));


    }
  }



  res = t(solve(info_Fish,score));
  return(res)
}

# % % Computation of E[d^2ln L/dbeta^2]: we use the formula \sum_{i,s,t}
# % % H_{i,s,t} X_{i,s} X_{i,t}', where H is the Hessian of denom
# %
# % for i=1:k
# %     for j=1:k
# %         for s=1:T
# %             XX(:,:,s)=bsxfun(@times,X,X(:,s));
# %             for t=1:T
# %                 if s~=t
# %                     temp = V(:,s).*V(:,t).*C_S_fun(S-2, V(:,(1:T)~=s & ...
#                                                             %                            (1:T)~=t))./denom;
# %                 else
#   %                     temp = J_denom(:,s);
#   %                 end
#   %                 hess_denom(:,s,t) = temp - J_denom(:,s).*J_denom(:,t);
#   %             end
#   %         end
#   %         hess(i,j) = - sum(XX.*hess_denom,[1 2 3])/n;
#   %     end
#   % end
#   %
#   % res = (hess \ score)';


