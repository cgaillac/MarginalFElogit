#' Approximation of the true effect and bounds on the ATE for the DGP in DDL
#'
#'
#' @param n_l the sample size used to compute the true bounds by Monte Carlo simulations (recommended 10e6)
#' @param case_DGP the number indexing the DGP in DDL
#' @param T the number of observed periods
#' @param beta0 a vector for the true values of the parameter beta
#' @param mat_C_k_T a matrix containing the combinatorial numbers
#' @param nbCores the number of cores used by the program to compute the true bounds on the AME. We strongly recommend to use nbCores>3.
#' @param selectX the selected covariate to evaluate the ATE
#' @param alpha_c a value to change the average of alpha. Default is 0.
#'
#' @return A list containing:
#'
#'  - the value of the true bounds on the AME;
#'
#'  - the true value of the AME.
#'
#' @export
#'
# @examples
true_param_bounds_paral_ATE <- function(n_l,case_DGP,T,beta0,mat_C_k_T,nbCores,selectX, alpha_c = 0 ){

    dimX= length(beta0)

    # Not necessary to simuulate if case_DGP=1
    if( case_DGP==1 & dimX ==1){

      true_param=  Lambda(beta0)-Lambda(0);
      true_bounds = c(true_param,true_param)
    }else{

      # n_l=1000000

      X= (array(rand(n_l,T*dimX),c(n_l,T,dimX))<=0.5)*1;
      # dim(  X)
      if( length(dim(X)) ==3){
        XT = matrix(X[,T,],n_l,dimX)
      }else{
        XT = matrix(X[,T],n_l,dimX)
      }

      Xbeta0 = zeros(n_l,T);

      # t= 1
      for( t in 1:T){
        Xbeta0[,t] = matrix(X[,t,],n_l,dimX)%*%beta0;
      }

      if(case_DGP==1){
        alpha = alpha_c+  zeros(n_l,1);
      }else if(case_DGP==2){
        alpha = alpha_c+  XT[,1] + 2*(rand(n_l,1)<=.5)-1;
      }else if(case_DGP==3){
        alpha = alpha_c+  XT[,1] + randn(n_l,1);
      }else{
        prob_un = pmax(0,abs(XT[,1])-1/4);
        temp = rand(n_l,1);
        alpha = alpha_c+  XT[,1] - (temp<=prob_un)*1 + (temp>1-prob_un)*1
      }

      index=Xbeta0+ alpha%*%rep(1,dim(X)[2])
      prob = Lambda(index);
      Y = (rand(n_l,T)<=prob)*1;

      V=exp(Xbeta0);

      # Matrix with k-th column = Pr(S=k|X=x)/C_k(x,b0)
      mat_ratio = zeros(n_l,T+1);

      ### change when dimX > 1
      XT0 = XT #matrix(0,n_l,dimX)
      XT0[XT[,selectX]==0,selectX ] = 1
      XT0[XT[,selectX]==1,selectX ] = 0


      term = Lambda(XT0%*%beta0+ alpha%*%rep(1,1))
      true_param = mean( Y[,T]*(2*XT[,selectX]-1)) + mean(term*(1-XT[,selectX]) - term*XT[,selectX] )

      # mean(XT==0)*mean(term[XT==0]) - mean(XT==1)*mean(term[XT==1])
      # true_param = mean( Y[,T]*(2*XT-1)) + mean(term*(XT==0) - term*(XT==1) )

      # % We first compute c_t (or an approximation of it) using
      # % c_t(x)=int_0^1 \frac{u^t}{\prod_{t=1}^{T-1} [1+u((exp((x_t-x_T)'
      #   % \beta_0) -1)]} dF_{U|X}(u|x) (see the proof of Lemma 1).

      Vtilde =  V/(matrix(exp(XT0%*%beta0),n_l,1)%*%rep(1,dim(V)[2]))

      # mat_p0 = .5*ones(n0,2)
      # # # % Probabilities associated with eta
      # case_DGP=2
      if( case_DGP==1){
        u = Lambda(XT0%*%beta0)
        mat_p =ones(n_l,1);

         }else if( case_DGP==2){

          u = cbind(Lambda(XT0%*%beta0+XT[,selectX]-1), Lambda(XT0%*%beta0+XT[,selectX]+1));
          mat_p = .5*ones(n_l,2);

        }else if (case_DGP==4){
          u = cbind(Lambda(XT0%*%beta0+XT[,selectX]-1), Lambda(XT0%*%beta0+XT[,selectX]),
                    Lambda(XT0%*%beta0+XT[,selectX]+1));

          P=pmax(abs(XT[,selectX])-1/4,0);
          mat_p=cbind(P,1-2*P,P);
        }else{
          # % DGP3: approximation of the integral by Gaussian quadrature
          nb_nodes=31;
          out = gausshermi(nb_nodes,0,1);
          w = out[[1]]
          roots= out[[2]]

          u=zeros(n_l,nb_nodes)
          mat_p=matrix(1,n_l,1)%*%t(w)
          temp=0;
          for(j in seq(1,nb_nodes)){
            u[,j] = Lambda(XT0%*%beta0+XT[,selectX]+roots[j])
          }
        }

      # % Computation of prod_{t=1}^{T-1} [1+u((exp((x_t-x_T)'\beta_0) -1)]
      # j=1
      denom= zeros(n_l,dim(u)[2]);
      for (j in 1:dim(u)[2]){
        Mat = u[,j]*(Vtilde-1)
        denom[,j] = apply(matrix(1+Mat,n_l,dim(Vtilde)[2]),1,prod)
      }

      c_true = zeros(n_l,T+1);
      # t=1
      for (t in 0:T){
        c_true[,t+1]= apply(mat_p*(u^t)/denom,1,sum);
      }
      #
      ###
      # T_minus_t = repmat(T - (0:T),n_l,1);
      # S = matrix(rowSums(Y))
      # C_S_vec = C_S_fun(S,Vtilde);
      # S_minus_t = matrix(S)%*%matrix(1,1,T+1) -   matrix(1,dim(matrix(S))[1],1)%*%(0:T);
      # mat_combin = choose(T_minus_t,S_minus_t)/repmat(C_S_vec,1,T+1);
      #
      # c_true = zeros(n_l,T+1);
      # rep = c(0,1)
      # xgrid = expand.grid(rep,rep)
      # # t=1
      #
      # for (t in 0:T){
      #   for(k in 1:dim(xgrid)[1]){
      #     c_true[rowSums(X[,,1]==(matrix(1,n_l,1)%*%as.matrix(xgrid[k,])))==T,t+1] = mean(mat_combin[rowSums(X[,,1]==(matrix(1,n_l,1)%*%as.matrix(xgrid[k,])))==T,t+1])
      #   }
      # }
      # c_true0 = c_true

      Vtilde_min1 = Vtilde - 1
      # Vtilde_min1_0 = cbind(V/exp(beta0) - 1)

      Mat_lambda=zeros(n_l,T+2);
      for (t in 0:T){
        Mat_lambda[,t+2]=C_S_fun(t,Vtilde_min1);
      }

      # Mat_lambda <- Mat_lambda*((1-XT)%*%rep(1,T+2)) + Mat_lambda1*(XT%*%rep(1,T+2))
      lambda_T_plus1 = matrix(Mat_lambda[,T+2], n_l,1)

      # % Identified term in Delta: expectation of the 1st term in Lemma 1
      Mat = matrix(Mat_lambda[,2:(T+1)]*c_true[,2:dim(c_true)[2]], n_l , T)

      term1= matrix(rowSums(Mat))
      term1 =  mean((1-XT[,selectX])*term1 - XT[,selectX]*term1)
      m_true = matrix(c_true[,2:dim(c_true)[2]]/c_true[,1], n_l, dim(c_true)[2]-1)

      # i =1
      nbCores=4
      sfInit(parallel = TRUE, cpus = nbCores, type = "SOCK")
      sfExportAll()
      sfLibrary(R.matlab)
      sfLibrary(pracma)
      res0 <- sfLapply(1:n_l, boot_MC,m_true,lambda_T_plus1,T,c_true)
      # res1 <- sfLapply(1:n1, boot_MC,m_true1,beta0 = matrix(1,dimX,1),lambda_T_plus1_1,T,c_true1)
      # boot_MC(400000,m_true,lambda_T_plus1,T,c_true)
      sfStop()
      # dim(lambda_T_plus1)
      # length(res0)

      # true_bounds1 = zeros(n1,2);
      true_bounds0 = zeros(n_l,2);
      for(i in 1:n_l){
        true_bounds0[i,1] <- res0[[i]][1]
        true_bounds0[i,2] <- res0[[i]][2]
      }


      # sum(is.na(m_true))
      # sum(is.na(true_bounds0[400000:600000,]))
      true_bounds =  mean( Y[,T]*(2*XT[,selectX]-1)) + term1 +
        mean(XT[,selectX]==0)* apply(true_bounds0[XT[,selectX]==0,],2,sumNA)/sum(XT[,selectX]==0) -
        mean(XT[,selectX]==1)* rev( apply(true_bounds0[XT[,selectX]==1,],2,sumNA)/sum(XT[,selectX]==1)) ;


      # true_bounds =  mean( Y[,T]*(2*XT-1)) +colMeans(  term2 *((1-XT)%*%rep(1,dim( term2 )[2])) -  term2 *(XT%*%rep(1,dim( term2 )[2]))) ;
      # true_param
      # true_bounds
      # colMeans(true_bounds*((XT==0)%*%rep(1,dim(true_bounds)[2])) - true_bounds*((XT==1)%*%rep(1,dim(true_bounds)[2])))
      # true_param = mean( Y[,T]*(2*XT-1)) + mean(XT==0)*  mean(term[XT==0]) - mean(XT==1)*  mean(term[XT==1]) ;
      # true_bounds =  mean( Y[,T]*(2*XT-1)) +  mean(XT==0)*( term1_0 +  colMeans(true_bounds[XT==0,])  ) - mean(XT==1)*  (term1_1 +c(B1[2],B1[1]) );
 }
    out=vector("list")
    out[[1]] <- true_bounds
    out[[2]] <- true_param
    return(out)
}
