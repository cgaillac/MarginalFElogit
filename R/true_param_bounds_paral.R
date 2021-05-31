#' Approximation of the true effect and bounds on the AME for the DGP in DDL
#'
#' @param n_l the sample size used to compute the true bounds by Monte Carlo simulations (recommended 10e6)
#' @param case_DGP the number indexing the DGP in DDL
#' @param T the number of observed periods
#' @param beta0 a vector for the true values of the parameter beta
#' @param mat_C_k_T a matrix containing the combinatorial numbers
#' @param nbCores the number of cores used by the program to compute the true bounds on the AME. We strongly recommend to use nbCores>3.
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
true_param_bounds_paral <- function(n_l,case_DGP,T,beta0,mat_C_k_T,nbCores , alpha_c = 0 ){

    dimX= length(beta0)

    # Not necessary to simuulate if case_DGP=1
    if( case_DGP==1 & dimX==1){
        true_param = 2*Lambda(.5)-1;
        true_bounds = cbind(true_param, true_param);
    }else{

      # X= array(randn(n_l,T*dimX),c(n_l,T,dimX));
      X= array(rand(n_l,T*dimX),c(n_l,T,dimX))-.5;
      # dim(  X)
      if( length(dim(X)) ==3){
        XT =  matrix(X[,T,],n_l,dimX)
      }else{
        XT = matrix(X[,T],n_l,dimX)
      }

      Xbeta0 = zeros(n_l,T);

      # t= 1
      for( t in 1:T){
        Xbeta0[,t] = matrix(X[,t,],n_l,dimX)%*%beta0;
      }

     # case_DGP=2
      if(case_DGP==1){
        alpha =  alpha_c+  zeros(n_l,1);
      }else if( case_DGP==2){
        alpha =   alpha_c+  XT[,1] + 2*(rand(n_l,1)<=.5)-1;
      }else if( case_DGP==3){
        alpha =  alpha_c+   XT[,1] + randn(n_l,1);
      }else{
        prob_un = pmax(0,abs(XT[,1])-1/4);
        temp = rand(n_l,1);
        alpha =   alpha_c+  XT[,1] - (temp<=prob_un) + (temp>1-prob_un);
      }

        # Necessary for case_DGP>1
        # We use a large sample of X's
      index=Xbeta0 + alpha%*%rep(1,dim(X)[2])
      V=exp(Xbeta0);

      # Matrix with k-th column = Pr(S=k|X=x)/C_k(x,b0)
      mat_ratio = zeros(n_l,T+1);


      true_param = matrix(beta0*mean(Lambda_prime(matrix(XT,n_l,dimX)%*%beta0+ alpha%*%rep(1,1)  )), dimX,1);

      # % We first compute c_t (or an approximation of it) using
      # % c_t(x)=int_0^1 \frac{u^t}{\prod_{t=1}^{T-1} [1+u((exp((x_t-x_T)'
      #   % \beta_0) -1)]} dF_{U|X}(u|x) (see the proof of Lemma 1).
      Vtilde =  V/(matrix(V[,T])%*%rep(1,dim(V)[2]))
      c_true = zeros(n_l,T+1);

      # u = cbind(Lambda(XT%*%beta0), Lambda(XT%*%beta0));
      # mat_p = .5*ones(n_l,2);
      if( case_DGP==1){
        u = cbind(Lambda(XT%*%beta0));
        # % Probabilities associated with eta
        mat_p = 1*ones(n_l,1);

      }else if( case_DGP==2){
          u = cbind(Lambda(XT%*%beta0+XT[,1]-1), Lambda(XT%*%beta0+XT[,1]+1));
          # % Probabilities associated with eta
          mat_p = .5*ones(n_l,2);

        }else if (case_DGP==4){
          u = cbind(Lambda(XT%*%beta0+XT[,1]-1), Lambda(XT%*%beta0+XT[,1]),Lambda(XT%*%beta0+XT[,1]+1));
          P=pmax(abs(XT[,1])-1/4,0);
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
            u[,j] = Lambda(XT%*%beta0+XT[,1]+roots[j])
          }
        }


      # % Computation of prod_{t=1}^{T-1} [1+u((exp((x_t-x_T)'\beta_0) -1)]
      # j=1
         if(T==1){
           denom= ones(n_l,dim(u)[2]);
         }else{
          denom= zeros(n_l,dim(u)[2]);
          for (j in 1:dim(u)[2]){
              Mat = u[,j]*(Vtilde[,1:(dim(Vtilde)[2]-1)]-1)
              denom[,j] = apply(matrix(1+Mat,n_l,(dim(Vtilde)[2]-1)),1,prod)
          }
         }
        # t=1
        for (t in 0:T){
          c_true[,t+1]= apply(mat_p*(u^t)/denom,1,sum);
        }

        Vtilde_min1 = cbind(Vtilde - 1, -ones(n_l,1))

        Mat_lambda=zeros(n_l,T+1);
        for (t in 1:T){
          Mat_lambda[,t+1]=C_S_fun(t-1,Vtilde_min1) ;
        }

        # ### could also be written
        #  Vtilde_min1 = Vtilde - 1
        # Mat_lambda0=zeros(n_l,T+1);
        # for (t in 1:T){
        #    Mat_lambda0[,t+1]=C_S_fun(t-1,Vtilde_min1[,1:(dim(Vtilde_min1)[2]-1)]) - C_S_fun(t-2,Vtilde_min1[,1:(dim(Vtilde_min1)[2]-1)]) ;
        # }

      # if((dim(Vtilde_min1)[2]-2)==1){
      #   lambda_T_plus1 = -apply(matrix(Vtilde_min1[,1:(dim(Vtilde_min1)[2]-2)], n_l,dim(Vtilde_min1)[2]-2),1,prod);
      # }else{
        if(T==1){
          lambda_T_plus1 = - matrix(1,n_l,1)
        }else{
          lambda_T_plus1 = -apply(matrix(Vtilde_min1[,1:(dim(Vtilde_min1)[2]-2)], n_l,dim(Vtilde_min1)[2]-2),1,prod);
         }
      # % Identified term in Delta: expectation of the 1st term in Lemma 1
      Mat = matrix(Mat_lambda[,2:dim(Mat_lambda)[2]]*c_true[,2:dim(c_true)[2]],n_l , dim(c_true)[2]-1)
      term1=beta0*sum(Mat)/(dim(Mat)[1]);


      m_true = c_true[,2:dim(c_true)[2]]/c_true[,1]

    # nbCores=4
    sfInit(parallel = TRUE, cpus = nbCores, type = "SOCK")
    sfExportAll()
    sfLibrary(R.matlab)
    sfLibrary(pracma)
    res0 <- sfLapply(1:n_l, boot_MC,m_true,lambda_T_plus1,T,c_true)
    sfStop()
    # boot_MC(1,m_true,lambda_T_plus1,T,c_true)

    true_bounds = zeros(n_l,2);
    for(i in 1:n_l){
        true_bounds[i,1] <- res0[[i]][1]
        true_bounds[i,2] <- res0[[i]][2]
    }
    # true_bounds  <- unlist()
    #
    #
    #
    true_bounds = matrix(term1,dimX,1)%*%rep(1,2)+ matrix(beta0,dimX,1)%*%colMeans(true_bounds)
 }
    out=vector("list")
    out[[1]] <- true_bounds
    out[[2]] <- true_param
    out[[3]] <- c_true
    return(out)
}
