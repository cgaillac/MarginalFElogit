#' Function which compute the ATE at a specified value of T contained in Tall.
#'
#' @param Yall matrix n x T containing the values of the dependent variable Yt.
#' @param Xall array of dimensions n x T x dimX containing the values of the predictors at the different periods.
#' @param prop_T a vector containin the proportions of the different numbers of periods observed in the dataset.
#' @param grid_T a vector containing the different numbers of periods observed in the dataset.
#' @param n_dist a vector containing the number of individuals which are observed for T periods, where T is given in grid_T
#' @param Tmax the maximal number of periods observed in the dataset
#' @param Tall vector of size n x1 containing the specified values of T where we compute the ATE.
#' @param Tinf vector of size n x1 containing the number of periods observed for each individual.
#' @param Call1 matrix n x 1 containing the identifiers of the clusters to which each individual belongs. Default is NULL
#' @param mat_C_k_T a matrix containing the combinatorial numbers.
#' @param cheb_coeff the coefficients of the Chebyshev polynomials T_n.
#' @param b_hat the estimated value of beta0 using the conditional maximum likelihood estimator.
#' @param alpha the confidence level for the confidence intervals. Default is 5\%.
#' @param CIOption the option for the choice of the type of confidence intervals for the quick method, either CI2 or CI3. Default is CI2.
#' @param Option  chosen method to compute the ATE: either sharp or quick. Default is quick.
#' @param dimX the number of covariates
#' @param selectX the selected covariate to compute the ATE (one by one)
#' @param phi_b the matrix of size n x dimX containing the inflence function of the  conditional maximum likelihood estimator of beta.
#' @param std_b the estimated standard deviation of the inflence function of the  conditional maximum likelihood estimator of beta.
#'
#' @return A list containing the values of : Option, n, Tmax, Time (computational time),
#'
#'  - Delta_hat: either the estimator of the bounds on Delta (sharp method) or the approximation of Delta (quick method);
#'
#'  - length_CI: the length of the confidence intervals;
#'
#'  - et: the estimated standard error of the influence function(s), either of the two bounds on Delta (sharp method) or of
#'  the approximation of Delta (quick method);
#'
#'  - bias_sup: the estimated upper bound on the bias (quick method);
#'
#'  - CI: the confidence intervals at the alpha level;
#'
#'  - phi_alpha: the phi alpha statistic intervening in the computation of the confidence interval of the sharp method;
#'
#'  - b_hat: the estimated value of beta0 using the conditional maximum likelihood estimator;
#'
#'  - std_b: the estimated standard deviation of the inflence function of the  conditional maximum likelihood estimator of beta;
#'
#'  - influence: the matrix of size n x dimX containing the inflence function of the  conditional maximum likelihood estimator of beta;
#'
#'  - for_average: a list containing the estimated values of Delta(x) and the influence function weighted by the values of Tinf, to compute the
#'  average of the AME over all the observed periods.
#' @export
#'
# @examples
compute_ATE_t <- function(Yall,Xall, prop_T, grid_T,n_dist,Tmax,Tall,Tinf,Call1,mat_C_k_T,cheb_coeff,b_hat,alpha, CIOption ,Option,dimX, selectX, phi_b, std_b){

  ### parameters for the sharp method.
  ratio=1
  RK = 1/(2*sqrt(pi)); #% int K^2(u)du
  kappa2 = 1; #% int u^2 K(u)du

  ### to stock the phi_alpha stat.
  phi_alpha_m= matrix(NA,dimX,1)

  ### number of variables to compute the AME/ATE
  if(is.null(selectX)){
    selectX = 1
  }

  start_time <- Sys.time()

  grid_T0 = sort(unique(Tall))
  nb_clust =length(unique(Call1))

  ## Estimation, either "quick" or "sharp"
  if(Option == "quick"){


    X = array(Xall[!is.na(Tall),,], c(sum(!is.na(Tall)),dim(Xall)[2],dim(Xall)[3]))
    Y = matrix(Yall[!is.na(Tall),], sum(!is.na(Tall)), dim(Yall)[2])
    Tall0 = Tall[!is.na(Tall)]
    Tinf0 = Tinf[!is.na(Tall)]
    Call10 = matrix(Call1[!is.na(Tall),], sum(!is.na(Tall)), 1)
    grid_T = grid_T[grid_T >=min( Tall0 )]
    n <- dim(X)[1]
    n_s <- n
    # T <- dim(X)[2]


    if(length(dim(X))==3){
      dimX = dim(X)[3]
    }else{
      dimX=1
    }

    S = apply(Y,1,sumNA)
    # dim(  X) ti=2
    XT =  matrix(NA,n,dimX)
    if( length(dim(X)) ==3){
      # XT = matrix(X[,Tall,],n,dimX)
      for (ti in grid_T0){
        if(sum(Tall0==ti)>0){
          XT[Tall0==ti,] =  matrix(X[Tall0==ti,ti,],sum(Tall0==ti),dimX)
        }
      }
    }else{
      # XT = matrix(X[,Tall],n,dimX)
      for (ti in grid_T0){
        if(sum(Tall0==ti)>0){
          XT[Tall0==ti,] = matrix(X[Tall0==ti,ti],sum(Tall0==ti),dimX)
        }
      }
    }

    if(length(dim(X))==2){
      Xtilde = X - matrix(XT,n,1)%*%rep(1,dim(X)[2])
    }else{
      Xtilde = X
      for(j in 1:dim(X)[3]){
        Xtilde[,,j] = X[,,j] - matrix(XT[,j],n,1)%*%rep(1,dim(X)[2])
      }
    }

    index = matrix(NA,n,Tmax)

    ## test if X of dim > 2
    if(length(dim(X))==3){
      for(ti in grid_T){
        for (t in 1:ti){
          index[Tinf0==ti,t] = X[Tinf0==ti,t,]%*%matrix(b_hat);
        }
      }
    }else{
      index = X*b_hat;
    }
    V = exp(index);

    # YT=Y[,T]*1;
    YT=matrix(NA, dim(Y)[1],1)
    for(ti in grid_T0){
      # YT = matrix(Y[Tall], dim(Y)[1],1)
      YT[Tall0==ti,] =Y[Tall0==ti,ti]
    }

    S = apply(Y,1,sumNA)

    ### first term in ATE
    term1=sum(YT*(2*XT[,selectX]-1))/n

    # % Coefficients of the Chebyshev polynomial
    # temp = coeff_Cheb01(T+1)/2^(T+1);
    # b = -rev(temp[2:length(temp)]);
    b = matrix(NA,n, Tmax+1)
    #ti=2
    for(ti in grid_T){
      temp = -cheb_coeff[[ti]][1:(ti+1)]/(2^(ti+1))
      b[Tinf0==ti,1:(ti+1)] <- matrix(1,sum(Tinf0==ti),1)%*%matrix(  temp, 1,ti+1)
      # deriv_a[Tall==ti,1:(ti+1)] = deriv_lambda[Tall==ti,1:(ti+1)] - matrix(deriv_lambda[Tall==ti,(ti+2)])%*%;
    }

    XT0 = XT #matrix(0,n_l,dimX)
    XT0[XT[,selectX]==0,selectX ] = 1
    XT0[XT[,selectX]==1,selectX ] = 0

    Vtilde =  V/(matrix(exp(XT0%*%b_hat),n,1)%*%rep(1,dim(V)[2]))
    Vtilde_T =  V/(matrix(exp(XT%*%b_hat),n,1)%*%rep(1,dim(V)[2]))

    # % Computation of the coefficients a_t(x, b_hat)
    mat_binom = matrix(NA,n,Tmax+1);
    lambda = matrix(NA,n,Tmax+2);
    lambda[,1] = 0
    Vmodif = Vtilde  - 1

    #ti=2
    # t=2
    for(ti in grid_T){
      for(t in 0:ti){
        mat_binom[Tinf0==ti,t+1] = choose((ti-t)*ones(sum(Tinf0==ti),1),S[Tinf0==ti]-t);
        lambda[Tinf0==ti,t+2] = C_S_fun(t,matrix(Vmodif[Tinf0==ti,1:ti], sum(Tinf0==ti), ti) );
      }
    }


    #### change for dimX >1
    f1 = matrix(NA,n,Tmax+1);
    f2 = matrix(NA,n,1);
    for(ti in grid_T){
      f1[Tinf0==ti,1:(ti+1)] = lambda[Tinf0==ti,1:(ti+1)]*(matrix(exp( b_hat[selectX]*S[Tinf0==ti])*(1-XT[Tinf0==ti,selectX]), sum(Tinf0==ti), 1))%*%rep(1,(ti+1))  -
        exp(- b_hat[selectX]*S[Tinf0==ti])*lambda[Tinf0==ti,1:(ti+1)]*(matrix(XT[Tinf0==ti,selectX], sum(Tinf0==ti),1)%*%rep(1,(ti+1)))
      f2[Tinf0==ti,1] = exp( b_hat[selectX]*S[Tinf0==ti])* lambda[Tinf0==ti,(ti+2)]*(1-XT[Tinf0==ti,selectX]) -
        exp(- b_hat[selectX]*S[Tinf0==ti])*lambda[Tinf0==ti,ti+2]*XT[Tinf0==ti,selectX]
    }

    ### matrix of coefficients d(x,s, b_hat)
    mat_a =  f1 + (f2%*%rep(1,Tmax+1)) * b;

    # Vtilde_T =  V/(matrix(exp(XT%*% b_hat),n,dimX)%*%rep(1,dim(V)[2]))
    vect_C_S  = matrix(NA,n,1)
    for(ti in grid_T){
      vect_C_S [Tinf0==ti] = C_S_fun(S[Tinf0==ti],matrix(Vtilde_T[Tinf0==ti,1:ti], sum(Tinf0==ti), ti));
    }



    ### change for dimX > 1.
    Mat_fin = mat_a*mat_binom/(vect_C_S%*%rep(1,dim(mat_a)[2]))
    # Mat_fin = mat_a*mat_binom*(exp(XT%*% b_hat*S)%*%rep(1,dim(mat_binom)[2]))/(vect_C_S%*%rep(1,dim(mat_a)[2]))
    term2 =sum(Mat_fin,na.rm=TRUE)/n; # mean(apply(Mat_fin, 1, sum  ));
    #
    Delta_hat = term1+term2;


    for_average = vector("list")
    for_average[[1]] =sum(YT*(2*XT[,selectX]-1)/Tinf0)/n+ sum(Mat_fin/(matrix(Tinf0,n,1)%*%rep(1,dim(mat_binom)[2])),na.rm=TRUE)/n;



    # true_param
    ### for the computation of the bias of  ATE_tilde
    for_bias = exp(- b_hat[selectX]*S) *abs(matrix(lambda[,Tall0],n,1))*XT[,selectX]+ exp( b_hat[selectX]*S) * abs(matrix(lambda[,Tall0],n,1))* (1-XT[,selectX]);
    bias_sup = mean(mat_binom[,1]*for_bias/vect_C_S/(2*4^Tinf0))

    for_average[[2]] =mean(mat_binom[,1]*for_bias/vect_C_S/(2*4^Tinf0*Tinf0))

    # phi = infl_func_beta(b_hat,Yall, Xall, Call1);
    # std_b =  apply(phi,2,std)/sqrt(dim(phi)[1])

    # Influence function of Delta_hat
    # Compute the influence function of beta_hat.
    # phi = infl_func_beta(b_hat,Y, X);


    if(dimX>1){
      term1 =  matrix(apply(Mat_fin,1,sumNA))+ YT*(2*XT[,selectX]-1)
    }else{
      term1 =  apply(Mat_fin,1,sumNA) + YT*(2*XT-1)
    }

    mat_combin = matrix(NA,n,max(grid_T)+1)
    for(ti in grid_T){
      T_minus_t = repmat(ti - (0:ti ),sum(Tinf0==ti),1);
      S_minus_t =  matrix(S[Tinf0==ti])%*%rep(1,ti +1) - matrix(rep(1,length(S[Tinf0==ti])))%*%seq(0,ti )  #bsxfun(@minus, S, seq(0,T));
      mat_combin[Tinf0==ti,1:(ti+1)] = choose(T_minus_t,S_minus_t)/repmat(matrix(vect_C_S[Tinf0==ti]),1,ti+1);
    }


    ## second term : has to compute derivative of 1) CS and 2) d from the one of lambda_t

    term2 = zeros(n,1);
    # k=1
    # ti=2
    # t=1
    for (k in (1:dimX)){
      # k= selectX
      deriv_C_S = matrix(0,n,1);
      if(length(dim(X)) == 2){
        for(ti in grid_T){
          for( t in 1:ti){

            Vtilde_T0 = matrix(Vtilde_T[Tinf0==ti,(1:ti)],sum(Tinf0==ti), ti )

            deriv_C_S[Tinf0==ti,] = deriv_C_S[Tinf0==ti,] +
              Xtilde [Tinf0==ti,t]*Vtilde_T[Tinf0==ti,t]*C_S_fun(S[Tinf0==ti]-1, matrix(Vtilde_T0[,(1:ti)!=t],sum(Tinf0==ti), sum((1:ti)!=t)  )  );
          }
        }
      }else{
        for(ti in grid_T){
          for( t in 1:ti){

            Vtilde_T0 = matrix(Vtilde_T[Tinf0==ti,(1:ti)],sum(Tinf0==ti), ti )

            deriv_C_S[Tinf0==ti,] = deriv_C_S[Tinf0==ti,] +
              Xtilde [Tinf0==ti,t,k]*Vtilde_T[Tinf0==ti,t]*C_S_fun(S[Tinf0==ti]-1, matrix(Vtilde_T0[,(1:ti)!=t],sum(Tinf0==ti), sum((1:ti)!=t)  )  );
          }
        }
      }

      deriv_e = matrix(0,n,Tmax+1);
      for(ti in grid_T){
        if((ti+2) <= (Tmax+1)){
          deriv_e[Tinf0==ti,(ti+2):(Tmax+1)] <- NA
        }
      }
      # deriv_e1 = zeros(n,T+1);

      if(length(dim(X))==2){
        Xtilde0 = X - matrix(XT0,n,1)%*%rep(1,dim(X)[2])
      }else{
        Xtilde0 = X
        for(j in 1:dim(X)[3]){
          Xtilde0[,,j] = X[,,j] - matrix(XT0[,j],n,1)%*%rep(1,dim(X)[2])
        }
      }


      for(ti in grid_T){
        Vmodif0 = Vmodif[,1:ti]
        if(length(dim( X)) == 2){
          for( t in 1:ti){
            for( u in 1:ti){

              Vmodif00 = matrix(Vmodif0[Tinf0==ti,1:ti ], sum(Tinf0==ti),ti)
              deriv_e[Tinf0==ti,t+1] = deriv_e[Tinf0==ti,t+1] + Xtilde0[Tinf0==ti,u]*Vtilde[Tinf0==ti,u]* C_S_fun(t-1, as.matrix(Vmodif00[,(1:ti)!=u ], dim(Vmodif00)[1], sum((1:ti)!=u)));
              # deriv_e1[,t+1] = deriv_e1[,t+1] + X[,u]*V[,u]* C_S_fun(t-1, as.matrix(Vmodif1[,(1:T)!=u]));
            }
          }
        }else{
          for( t in 1:ti){
            for( u in 1:ti){

              Vmodif00 = matrix(Vmodif0[Tinf0==ti,1:ti ], sum(Tinf0==ti),ti)
              deriv_e[Tinf0==ti,t+1] = deriv_e[Tinf0==ti,t+1] + Xtilde0[Tinf0==ti,u,k]*Vtilde[Tinf0==ti,u]* C_S_fun(t-1, matrix(Vmodif00[,(1:ti)!=u], dim(Vmodif00)[1], sum((1:ti)!=u)));
              # deriv_e1[,t+1] = deriv_e1[,t+1] + X[,u,k]*V[,u]* C_S_fun(t-1, as.matrix(Vmodif1[,(1:T)!=u]));
            }
          }
        }
      }


      deriv_lambda = matrix(NA,n,Tmax+2);
      # deriv_lambda1 = zeros(n,T+2);
      for(ti in grid_T){
        deriv_lambda[Tinf0==ti,1] = 0
        deriv_lambda[Tinf0==ti,2:(ti+2)] = deriv_e[Tinf0==ti,1:(ti+1)]
      }
      # deriv_lambda1[,2:(T+2)] = deriv_e1[,1:(T+1)]


      if(k== selectX){

        f1 = matrix(NA,n,Tmax+1)
        f2 = matrix(NA,n,1)
        for(ti in grid_T){
          # Derivative of the coeff a_t wrt beta
          f1[Tall0==ti,1:(ti+1)]= (matrix(exp( b_hat[selectX]*S[Tall0==ti])*(1-XT[Tall0==ti,selectX]))%*%rep(1,(ti+1)))*(deriv_lambda[Tall0==ti,1:(ti+1)] +
                                                (matrix(S[Tall0==ti],sum(Tall0==ti),1)%*%rep(1,(ti+1)))*lambda[Tall0==ti,1:(ti+1)]) -
            (matrix(exp(- b_hat[selectX]*S[Tall0==ti])*XT[Tall0==ti,selectX])%*%rep(1,(ti+1)))*(deriv_lambda[Tall0==ti,1:(ti+1)] -
                                            (matrix(S[Tall0==ti],sum(Tall0==ti),1)%*%rep(1,(ti+1)))*lambda[Tall0==ti,1:(ti+1)])
          # f1p = lambda[,1:(T+1)]*(S*(exp( b_hat*S)*(1-XT))%*%rep(1,(T+1)))
          # f2 = deriv_lambda1[,1:(T+1)]*(XT%*%rep(1,(T+1)))
          f2[Tall0==ti,1] = exp( b_hat[selectX]*S[Tall0==ti])* (deriv_lambda[Tall0==ti,(ti+2)] + S[Tall0==ti]*lambda[Tall0==ti,(ti+2)])*(1-XT[Tall0==ti,selectX]) -
            exp(- b_hat[selectX]*S[Tall0==ti])*XT[Tall0==ti,selectX]*(deriv_lambda[Tall0==ti,ti+2]- S[Tall0==ti]*lambda[Tall0==ti,ti+2])
        }


      }else{

        f1 = matrix(NA,n,Tmax+1)
        f2 = matrix(NA,n,1)
        for(ti in grid_T){
          # Derivative of the coeff a_t wrt beta
          f1[Tall0==ti,1:(ti+1)]= (matrix(exp( b_hat[selectX]*S[Tall0==ti])*(1-XT[Tall0==ti,selectX]))%*%rep(1,(ti+1)))*deriv_lambda[Tall0==ti,1:(ti+1)]  -
            (matrix(exp(- b_hat[selectX]*S[Tall0==ti])*XT[Tall0==ti,selectX])%*%rep(1,(ti+1)))*deriv_lambda[Tall0==ti,1:(ti+1)]

          f2[Tall0==ti,1] = exp( b_hat[selectX]*S[Tall0==ti])* (deriv_lambda[Tall0==ti,(ti+2)])*(1-XT[Tall0==ti,selectX]) -
            exp(- b_hat[selectX]*S[Tall0==ti])*XT[Tall0==ti,selectX]*(deriv_lambda[Tall0==ti,ti+2])
        }


      }
      ### matrix of coefficients d(x,s, b_hat)
      deriv_a =  f1  + (f2%*%rep(1,dim(b)[2])) * b;
      # deriv_a = deriv_lambda[,1:(T+1)] - matrix(deriv_lambda[,(T+2)])%*%cheb_coeff[1:(T+1)]/(2^(T+1));

      ### add the other term for dimX > 1
      # mat_derivF = -  mat_a* ((deriv_C_S/vect_C_S^2)%*%rep(1,dim(mat_a)[2])) +  #mat_a*((exp( b_hat*S)/vect_C_S* b_hat*S)%*%rep(1,dim(mat_a)[2])) +
      #   ((1/vect_C_S)%*%rep(1,dim(mat_a)[2]))* deriv_a
      # # moy_deriv = (T+1)*mean(mat_combin*(deriv_a - bsxfun(@times, Mat_a, deriv_C_S/C_S_vec)),'all');
      # moy_deriv = mean(apply(mat_binom*  mat_derivF, 1, sum));
      #

      moy_deriv = sum(mat_combin*(deriv_a -  mat_a* ((deriv_C_S/vect_C_S)%*%rep(1,dim(mat_a)[2]))),na.rm=T)/n;
      term2 = term2 + moy_deriv*phi_b[!is.na(Tall),k];

      # term2 = b_hat * moy_deriv * phi;

    } ## end of the for k loop


    infl =  term1 + term2
    # et = std(infl)/sqrt(dim(infl)[1])

    nb_clust =length(unique(Call10))
    civals = unique(Call10)

    et = matrix(0,1,dim(infl)[2])
    for(ci0 in 1:nb_clust){
      # if(dimX==1){
        ci = civals[ci0]
      # if(dimX==1){
        et =et + var(infl[Call10==ci ,])
      # }else{
        # et =et + diag(var(infl[Call10==ci ,]))
      # }
    }
    et =sqrt(et)/sqrt(dim(infl)[1])

    eps_n = (2*log(log(n)))^(1/2)/sqrt(dim(infl)[1])

    if( CIOption == "CI3"){
      # length_CI[s] = 2 * et[s] * sqrt(ncx2inv(0.95,1,(bias_sup[s]/et[s])^2));
      length_CI = 2 * et * sqrt(qchisq(0.95,1, ncp = ((bias_sup + eps_n)/et)^2, lower.tail = TRUE, log.p = FALSE))
    }else{
      length_CI = 2 * et * sqrt(qchisq(0.95,1, ncp = (bias_sup/et)^2, lower.tail = TRUE, log.p = FALSE))
    }
    ### compute CI

    CI =  matrix(Delta_hat, 1,1)%*%rep(1,2)
    CI[,1] = CI[,1] -  matrix(length_CI/2, 1,1)
    CI[,2] = CI[,2] +  matrix(length_CI/2, 1,1)






  }else{ ### end of "option" quick



    nb_clust =length(unique(Call1))
    civals = unique(Call1)

    n <- dim(Yall)[1]
    n_s <- n
    # "option" slow
    boundsall = vector("list")
    infl = matrix(0,n,2)


    # phi_b =matrix(0,n,dimX)


    Tall0 = Tall[!is.na(Tall)]
    grid_T = grid_T[grid_T >=min( Tall0 )]
    n_s <- sum(!is.na(Tall))


    for_average = vector("list")

    # ci =1
    # ti0=1
    count = 1
    for(ci0 in 1:nb_clust){
      for(ti0 in 1:length(grid_T)){


        T = grid_T[ti0]
        ci = civals[ci0]
        selecti = Tinf==T & Call1==ci & !is.na(Tall)

        Y = matrix(Yall[selecti,1:T], sum(selecti) , T )
        X = array(Xall[selecti, 1:T ,], c(sum(selecti), T ,dimX))
        Tall0 <- Tall[selecti]
        n <- dim(X)[1]
        T <- dim(X)[2]

        # % Compute the influence function of beta_hat. Useful for inference on
        # % Delta, at the end.
        phi  = phi_b[selecti,]


        if(length(dim(X))==3){
          dimX = dim(X)[3]
        }else{
          dimX=1
        }

        T_minus_t = repmat(T - (0:T),n,1);
        S = apply(Y,1,sum)

        # YT=Y[,T]*1;
        YT=matrix(NA, dim(Y)[1],1)
        for(ti in grid_T0){
          if(sum(Tall0==ti)>0){
            # YT = matrix(Y[Tall], dim(Y)[1],1)
            YT[Tall0==ti,] =Y[Tall0==ti,ti]
          }
        }

        # dim(  X)
        # if( length(dim(X)) ==3){
        #   XT =  matrix(X[,T,],n,dimX)
        # }else{
        #   XT = matrix(X[,T],n,dimX)
        # }
        XT =  matrix(NA,n,dimX)
        if( length(dim(X)) ==3){
          # XT = matrix(X[,Tall,],n,dimX)
          for (ti in grid_T0){
            if(sum(Tall0==ti)>0){
              XT[Tall0==ti,] =  matrix(X[Tall0==ti,ti,],sum(Tall0==ti),dimX)
            }
          }
        }else{
          # XT = matrix(X[,Tall],n,dimX)
          for (ti in grid_T0){
            if(sum(Tall0==ti)>0){
              XT[Tall0==ti,] = matrix(X[Tall0==ti,ti],sum(Tall0==ti),dimX)
            }
          }
        }

        # ### change when dimX > 1
        # XT0 = matrix(0,n,dimX)
        # XT0[XT==0] = 1

        XT0 = XT #matrix(0,n_l,dimX)
        XT0[XT[,selectX]==0,selectX ] = 1
        XT0[XT[,selectX]==1,selectX ] = 0


        # Vtilde =  V/(matrix(exp(XT0%*% b_hat),n,dimX)%*%rep(1,dim(V)[2]))


        # %% Simulations
        # % Initialisation

        # std_bounds = zeros(1,2)
        # cov_rate = zeros(1,1)
        #
        # # % Percent of estimated I(m) equal to the true I(m)
        # tx_class=zeros(1,1);
        #
        # # % Parameter used to compute numerical derivatives wrt beta and gamma
        hderiv=1e-3;
        #

        # % Initialisation
        C_mat = zeros(n,T+1);
        PSt_X = zeros(n,T+1);
        # % Step 1: CMLE
        # options=optimoptions('fminunc','Display','off');
        # b_hat = optim(par =  b_hat, log_lik_FE_logit,Y=Y,X=X)$par


        # % Compute the influence function of beta_hat. Useful for inference on
        # % Delta, at the end.
        # phi = infl_func_beta(b_hat,Y, X);
        # phi_b[Tinf==T & Call1==ci,] =  phi


        index = zeros(n,T)

        ## test if X of dim > 2
        if(length(dim(X))==3){
          for ( t in 1:T){
            index[,t] = X[,t,]%*%matrix(b_hat);
          }
        }else{
          index = X*b_hat;
        }


        # % Intercept assuming that Y_T on X_T is a logit: used for the
        # % computation of the bandwidth below
        if(length(dim(X))==3){
          intercept = estim_intercept(Y,index);
        }else{
          intercept = estim_intercept(Y,index);
        }

        index_int = intercept + index;
        V = exp(index);

        # % Step 2: bounds on the AME

        # % We first need nonparametric estimators of mt(x).
        # % We use mt(x)=ct(x)/c0(x) for t\geq 0, with (c0(x),...,cT(x))'
        # % = A x (Pr(S=0|X=x)/C0(x; b_hat),...,Pr(S=T|X=x)/CT(x; b_hat).



        # %% Computation of the bandwidth for local linear estimation.
        # % We fix h so that estimated avg absolute bias = ratio * estimated
        # % avg asymp. std dev. To estimate both, we assume that the fixed
        # % effect is actually constant (and then estimated by intercept).


        ####################################   SELECT according to support and dimension of X
        # % Uniform draws on Supp(X). Used below to get the bandwidth.

        n_mx = 3000
        mx= array(0,c(n_mx,T,dimX))
        scaling = matrix(1,dimX,T)
        for(k in 1:dimX){
          # Fx = matrix(quantile(X[,,k],  probs = rand( n_mx,T)) ,n_mx,T)
          # mx[,,k] = Fx
          scaling[k,] = (apply(matrix(X[,,k],n,T),2,max)-apply(matrix(X[,,k],n,T),2,min))
          mx[,,k] = as.matrix(rep(1,n_mx))%*%apply(matrix(X[,,k],n ,T ),2,min)+ (as.matrix(rep(1,n_mx))%*%scaling[k,]) * rand(n_mx,T)
        }

        scaling0 = prod(c(scaling))

        randunif_x =mx;
        index_bis= matrix(0,dim(randunif_x[,,1])[1],dim(randunif_x[,,1])[2])
        for(k in 1:dimX){
          index_bis = index_bis + randunif_x[,,k]  * b_hat[k];
          # index_bis = index_bis + (randunif_x[,,k] >0.5) * b_hat[k];
        }

        # t = 2
        h = zeros(1,T+1);
        for (t in 0:T){
          C_mat[,t+1] = C_S_fun(t,V);
          # % Computation of the avg std dev: (R(k)^q int sigma^2(x)dx/n)^{1/2}
          # % (see Hansen's notes). We use the uniform draws above. The
          # % estimator is an average of sigma^2(x) over these draws
          mean_approx_St = exp(t*intercept)*apply(1-Lambda(intercept+index_bis),1,prod)* C_S_fun(t,exp(index_bis));
          std_approx_St = scaling0 *sqrt(RK^(dimX*T) * mean(mean_approx_St*(1-mean_approx_St))/n);

          # % Computation of the average absolute bias
          mean_approx_St = exp(t*intercept)*apply(Lambda(index_int),1,prod)*C_mat[,t+1];

          d0c = dim(Lambda(index_int))[2]
          d0r = dim(Lambda(index_int))[1]

          term1 = matrix(0,d0r, d0c)
          term2 = matrix(0,d0r, d0c)


          if (t==0){
            deriv_log_C = 0;
            for(k in 1:dimX){
              term1 = term1 +  (deriv_log_C + b_hat[k] * Lambda(index_int))^2;
              term2 = term2 + deriv_log_C*(b_hat[k] - deriv_log_C) + b_hat[k]^2*Lambda(index_int,1);
            }
          }else{
            bbx = V* (matrix(C_mat[,t])%*%rep(1,dim(V)[2])) #bsxfun(@times, V,  C_mat[,t])
            for(k in 1:dimX){
              deriv_log_C = b_hat[k] * ( bbx / (matrix(C_mat[,t+1])%*%rep(1,dim(bbx)[2])))  # bsxfun(@rdivide, ,C_mat[,t+1]);
              term1 =  term1 + (deriv_log_C + b_hat[k] * Lambda(index_int))^2;
              term2 =  term2 + deriv_log_C*(b_hat[k] - deriv_log_C) + b_hat[k]^2*Lambda(index_int,1);
            }
          }

          bias_approx_St = sum(mean(abs( (as.matrix(mean_approx_St)%*%rep(1,dim(term1)[2]))*(term1+term2))));

          # % Choice of the bandwidth
          h[1,t+1] = (std_approx_St / (sqrt(ratio)*bias_approx_St))^(1/(2+dimX*T/2));
        }
        # h= h/2

        # t=2
        for (t in 0:T){
          St = matrix((S==t)*1);
          if(t==T){
            out =  local_lin(St,X,h[,t+1]);
            ES_T = out[[1]]
            f= out[[2]]
          }else{
            out = local_lin(St,X,h[,t+1]);
            ES_T = out[[1]]
          }
          PSt_X[,t+1] = pmin(pmax(ES_T,0),1);
        }


        # % Normalization: to be sure that the probabilities sum to 1.
        # PSt_X = bsxfun(@rdivide, PSt_X, sum(PSt_X,2));
        PSt_X = PSt_X / (matrix(apply(PSt_X,1,sum))%*%rep(1,dim(PSt_X)[2]));

        T_minus_t = repmat(T - (0:T),n,1);
        deriv_beta = zeros(dimX,2);
        deriv_gamma = zeros(T+1,2);
        # version ="new"
        # s=1
        # if(version == "new"){
        #
        # # % Computation of the bounds and Pr(Ihat(m_i)=I(m_i))
        out =  Delta_bounds_ATE(b_hat, V, S, dimX, PSt_X, f, h,mat_C_k_T[[T]], RK, T_minus_t, XT0, XT, YT,selectX,Tall0)
        boundsall[[count]] = out[[1]]
        bounds = out[[1]]
        bounds_ind = out[[2]]

        # % Computation of the influence function



        for (i in 1:dimX){
          out = Delta_bounds_ATE(b_hat+hderiv*((1:dimX)==i),V, S, dimX, PSt_X, f, h,mat_C_k_T[[T]], RK, T_minus_t, XT0, XT, YT,selectX,Tall0);
          for_deriv_b= out[[1]]
          # bounds_ind = out[[2]]
          # tx_class[s] = out[[3]]
          deriv_beta[i,] = (for_deriv_b - bounds)/hderiv;
        }


        for (t in 0:T){
          out = Delta_bounds_ATE(b_hat, V, S, dimX, PSt_X+hderiv*ones(n,1)%*%((0:T)==t), f, h, mat_C_k_T[[T]], RK, T_minus_t, XT0, XT, YT,selectX,Tall0);
          for_deriv_g <- out[[1]]
          deriv_gamma[t+1,] = (for_deriv_g - bounds)/hderiv;
        }


        # }else{

        # % Computation of the bounds and Pr(Ihat(m_i)=I(m_i))
        # out = Delta_bounds_old(b_hat, X, PSt_X, f, h, mat_C_k_T, RK, case_DGP);
        # bounds[s,] = out[[1]]
        # bounds_ind = out[[2]]
        # tx_class[s] = out[[3]]
        #
        # temps[s,1]= start -stop;
        #
        #
        #
        # # % Computation of the influence function
        #
        # for (i in 1:dimX){
        #   out = Delta_bounds_old(b_hat+hderiv*((1:dimX)==i), X, PSt_X, f,  h, mat_C_k_T, RK, case_DGP);
        #   for_deriv_b= out[[1]]
        #   # bounds_ind = out[[2]]
        #   # tx_class[s] = out[[3]]
        #   deriv_beta[i,] = (for_deriv_b - bounds[s,])/hderiv;
        # }
        #
        #
        # for (t in 0:T){
        #   out = Delta_bounds_old(b_hat, X, PSt_X+hderiv*ones(n,1)%*%((0:T)==t),f, h, mat_C_k_T, RK, case_DGP);
        #   for_deriv_g <- out[[1]]
        #   deriv_gamma[t+1,] = (for_deriv_g - bounds[s,])/hderiv;
        # }
        # } ## end of the if "version"


        # % Term corresponding to the influence of the nonparametric
        # % estimation of the Pr(S=t|X)
        Z_EZX = matrix(S)%*%rep(1,T+1)==(matrix(rep(1,length(S)))%*%t(matrix(0:T))) -  PSt_X;

        # infl = rbind(infl, bounds_ind + phi %*% deriv_beta + Z_EZX %*% deriv_gamma)
        infl[selecti,]  = bounds_ind + phi %*% deriv_beta + Z_EZX %*% deriv_gamma
        count=count+1


        sfStop()


      } ### end of T_i loop
    }

    count=1
    Delta_hat= matrix(0,1,2)
    for(ci0 in 1:nb_clust){
      for(ti in grid_T){

        ci = civals[ci0]
        selecti = Tinf==T & Call1==ci & !is.na(Tall)

        Delta_hat=  Delta_hat+ boundsall[[count]]*sum(selecti)/n_s
        count=count+1
      }
    }

    for_average = vector("list")
    count=1
    Delta_hat0= matrix(0, 1,2)
    for(ci0 in 1:nb_clust){
      for(ti in grid_T){

        ci = civals[ci0]
        selecti = Tinf==T & Call1==ci & !is.na(Tall)

        Delta_hat0=  Delta_hat0 + boundsall[[count]]*sum(selecti)/(n_s*ti)
        count=count+1
      }
    }
    for_average[[1]] <- Delta_hat0


    # Delta_hat
    et =apply(infl[!is.na(Tall) ,],2,std)/sqrt(n_s);
    # length_CI = Delta_hat[2] +et[2] - Delta_hat[1] +et[1]
    bias_sup = matrix(0,2,1)

    # stop =  Sys.time()
    # time =  start_time -stop;

    CI = matrix(NA,1,2)

    ## compute CI at level alpha.
    # k=1
    # for(k in 1:dimX){
    quant = quantile_IM(alpha, Delta_hat, et);

    # phi_alpha=1{sqrt(n)|beta hat_k|/std(phi_k)>q_{1-alpha/2}}
    phi_alpha = sqrt(n_s)*abs(b_hat[selectX])/std(infl[!is.na(Tall) ,])> qnorm(1-alpha/2)
    if(phi_alpha){
      CI = cbind(Delta_hat[1] - quant* et[1], Delta_hat[2] + quant*et[2]);
    }else{
      CI = cbind(min(0,Delta_hat[1] - quant* et[1]), max(0,Delta_hat[2] + quant*et[2]));
    }

    phi_alpha_m = phi_alpha*1
    length_CI = CI[,2] - CI[,1]
    # }
    std_b = apply(phi_b,2,std)/sqrt(dim(phi_b)[1])
  }

  # Save results
  end_time <- Sys.time()
  time = end_time - start_time


  out <- vector("list")
  out[[1]] <- Option
  out[[2]] <-  n_s
  out[[3]] <- grid_T
  out[[4]] <- as.numeric(time)
  out[[5]] <- Delta_hat
  out[[6]] <- length_CI
  out[[7]] <- et
  out[[8]] <- bias_sup
  out[[9]] <- CI
  out[[10]] <- phi_alpha_m
  out[[11]] <- b_hat
  out[[12]] <- std_b
  out[[13]] <- infl
  out[[14]] <-  for_average

  names(out) <- c("Option","n","Tmax","Time","Delta_hat","length_CI","et","bias_sup","CI","phi_alpha","b_hat","std_b", "influence","for_average")



  # computation of the influence function. + ASN.

  return(out)
}
