#' Function which compute the AME at a specified value of T contained in Tall.
#'
#' @param Yall matrix n x T containing the values of the dependent variable Yt.
#' @param Xall array of dimensions n x T x dimX containing the values of the predictors at the different periods.
#' @param prop_T a vector containin the proportions of the different numbers of periods observed in the dataset.
#' @param grid_T a vector containing the different numbers of periods observed in the dataset.
#' @param n_dist a vector containing the number of individuals which are observed for T periods, where T is given in grid_T
#' @param Tmax the maximal number of periods observed in the dataset
#' @param Tall vector of size n x1 containing the specified values of T where we compute the AME.
#' @param Tinf vector of size n x1 containing the number of periods observed for each individual.
#' @param Call1 matrix n x 1 containing the identifiers of the clusters to which each individual belongs. Default is NULL
#' @param mat_C_k_T a matrix containing the combinatorial numbers.
#' @param cheb_coeff the coefficients of the Chebyshev polynomials T_n.
#' @param b_hat the estimated value of beta0 using the conditional maximum likelihood estimator.
#' @param alpha the confidence level for the confidence intervals. Default is 5\%.
#' @param CIOption the option for the choice of the type of confidence intervals for the quick method, either CI2 or CI3. Default is CI2.
#' @param Option  chosen method to compute the AME: either sharp or quick. Default is quick.
#' @param dimX the number of covariates
#' @param selectX the vector of selected covariates to compute the AME.
#' If NULL then bounds are computed for all covariates. NULL by default.
#' @param phi_b the matrix of size n x dimX containing the influence function of the  conditional maximum likelihood estimator of beta.
#' @param std_b the estimated standard deviation of the influence function of the  conditional maximum likelihood estimator of beta.
#' @param nbCores the number of cores used by the program to compute the AME for the ``sharp" method.
#' @param ratio the ratio R in DDL for the nonparametric estimator of the conditional moments of S
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
#'  - std_b: the estimated standard deviation of the influence function of the  conditional maximum likelihood estimator of beta;
#'
#'  - influence: the matrix of size n x dimX containing the influence function of the  conditional maximum likelihood estimator of beta;
#'
#'  - for_average: a list containing the estimated values of Delta(x) and the influence function weighted by the values of Tinf, to compute the
#'  average of the AME over all the observed periods.
#' @export
#'
# @examples
compute_AME_t<- function(Yall,Xall, prop_T, grid_T,n_dist,Tmax,Tall,Tinf,Call1,mat_C_k_T,cheb_coeff,b_hat,alpha = 0.05, CIOption ,Option,dimX,selectX=NULL,  phi_b, std_b, nbCores=4,  ratio=10){

  ### parameters for the sharp method.
  # ratio=10
  RK = 1/(2*sqrt(pi)); #% int K^2(u)du
  kappa2 = 1; #% int u^2 K(u)du

  c_hat = vector("list")

  ### number of variables to compute the AME/ATE
  if(is.null(selectX)){
    nb_var = dimX
    selectX = 1:dimX
  }else{
    nb_var = length(selectX)
  }

  ### to stock the phi_alpha stat.
  phi_alpha_m= matrix(NA,nb_var,1)

  ###
  grid_T0 = sort(unique(Tall))
  start_time <- Sys.time()

  # dim(Xall)
  ## Estimation, either "quick" or "sharp"
  if(Option == "quick"){

    X = array(Xall[!is.na(Tall),,], c(sum(!is.na(Tall)),dim(Xall)[2],dim(Xall)[3]))
    Y = matrix(Yall[!is.na(Tall),], sum(!is.na(Tall)), dim(Yall)[2])
    Tall0 = Tall[!is.na(Tall)]
    Tinf0 = Tinf[!is.na(Tall)]
    Call10 = matrix(Call1[!is.na(Tall),], sum(!is.na(Tall)), 1)
    grid_T1 = grid_T[grid_T >=min( Tall0 )]

    # Tmax = max(grid_T1)
    # X = Xall
    # Y = Yall
    n <- dim(X)[1]
    # T <- dim(X)[2]
    ## stock the sample size
    n_s <- n

    if(length(dim(X))==3){
      dimX = dim(X)[3]
    }else{
      dimX=1
    }

    S = apply(Y,1,sumNA)
    # dim(  X)
    # if( length(dim(X)) ==3){
    #   XT =  matrix(X[,Tall,],n,dimX)
    # }else{
    #   XT = matrix(X[,Tall],n,dimX)
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

    # Step 1: CMLE
    # b_hat = fminunc(@(b) log_lik_FE_logit(b,Y,X), beta0,options);
    # if(dim(X)[2]==1){
    #   b_hat = optimize( log_lik_FE_logit,start_bounds,Y=Y,X=X)$minimum
    # }else{
    #   b_hat = optim(par = start_point, log_lik_FE_logit,Y=Y,X=X)$par
    # }
    # b_hat$par


    index = matrix(NA,n,Tmax)

    ## test if X of dim > 2
    if(length(dim(X))==3){
      for(ti in grid_T1){
        for (t in 1:ti){
          index[Tinf0==ti,t] = X[Tinf0==ti,t,]%*%matrix(b_hat);
        }
      }
    }else{
      index = X*b_hat;
    }
    V = exp(index);


    Vtilde = matrix(NA,dim(V)[1],dim(V)[2])
    for(ti in grid_T1){
      Vtilde[Tinf0==ti,] =  V[Tinf0==ti,]/(matrix(V[Tinf0==ti,Tall0[Tinf0==ti]],sum(Tinf0==ti),1)%*%rep(1,dim(V)[2])) #bsxfun(@rdivide,V,V[,T]);
    }

    Vtilde_min1 =  matrix(NA,dim(V)[1], (Tmax -1) )

    # ti=4
    for(ti in grid_T1){
      if(ti>=2){
        ind = 1:ti
        Vtilde_min1[Tinf0==ti,1:(ti-1)] =  Vtilde[Tinf0==ti, ind[ind!=Tall0[Tinf0==ti][1]]] - 1;
      }
    }


    out=best_approx_poly(Vtilde_min1,grid_T1,Tinf0);
    res = out[[1]]
    g= out[[2]]


    if(length(dim(X))==2){
      Xtilde = X - matrix(XT,n,1)%*%rep(1,dim(X)[2])
    }else{
      Xtilde = X
      for(j in 1:dim(X)[3]){
        Xtilde[,,j] = X[,,j] - matrix(XT[,j],n,1)%*%rep(1,dim(X)[2])
      }
    }



    # ti=3
    C_S_vec = matrix(NA,n,1)
    for(ti in grid_T1){
      C_S_vec[Tinf0==ti] = C_S_fun(S[Tinf0==ti],matrix(Vtilde[Tinf0==ti,1:ti], sum(Tinf0==ti) ,ti ));
    }


    mat_combin = matrix(NA,n,max(grid_T1)+1)
    for(ti in grid_T1){
      T_minus_t = repmat(ti - (0:ti ),sum(Tinf0==ti),1);
      S_minus_t =  matrix(S[Tinf0==ti])%*%rep(1,ti +1) - matrix(rep(1,length(S[Tinf0==ti])))%*%seq(0,ti )  #bsxfun(@minus, S, seq(0,T));
      mat_combin[Tinf0==ti,1:(ti+1)] = choose(T_minus_t,S_minus_t)/repmat(matrix(C_S_vec[Tinf0==ti]),1,ti+1);
    }

    Mat_fin  = matrix(NA,n,max(grid_T1)+1)
    Mat_a  = matrix(NA,n,max(grid_T1)+1)
    # Mat_fin = Mat_a*mat_combin;
    # ti =3
    for(ti in grid_T1){
      Mat_a[Tinf0==ti,1:(ti+1)] = fliplr(matrix(res[Tinf0==ti,1:(ti+1)], sum(Tinf0==ti),ti+1));
      Mat_fin[Tinf0==ti,1:(ti+1)] = Mat_a[Tinf0==ti,1:(ti+1)] *mat_combin[Tinf0==ti,1:(ti+1)];
    }


    moy_fin =  sum(Mat_fin,na.rm=T)/n;
    Delta_hat = b_hat[selectX] * moy_fin;

    for_average = vector("list")
    for_average[[1]] =  b_hat[selectX] * sum((Mat_fin)/(matrix(Tinf0,n,1)%*%rep(1,dim(mat_combin)[2])),na.rm=T)/n;


    ###############################
    m0 = mat_combin[,1]
    m0[is.na(m0)] = 0
    g[is.na(g)] = 0

    # Estimation of the maximal bias
    bias_sup = abs(b_hat[selectX]) * c(t(m0)%*%abs(g/4^Tinf0)/(n*2));
    for_average[[2]] = abs(b_hat[selectX]) * c(t( m0)%*%abs(g/(4^Tinf0*Tinf0))/(n*2));

    # Influence function of Delta_hat
    if(dimX>1){
      term1 =  matrix(apply(Mat_fin,1,sumNA))%*%b_hat + moy_fin *   phi_b[!is.na(Tall),]
    }else{
      term1 =  apply(Mat_fin,1,sumNA)*b_hat + moy_fin *   phi_b[!is.na(Tall),]
    }

    for_term2 = zeros(n,1);

    # k=1
    for (k in (1:dimX)){
      deriv_C_S = zeros(n,1);
      if(length(dim( Xtilde)) == 2){
        for(ti in grid_T1){
          if(ti==1){
            t=1
            deriv_C_S[Tinf0==ti] = deriv_C_S[Tinf0==ti] + Xtilde[Tinf0==ti,t]*Vtilde[Tinf0==ti,t]
          }else{
            for( t in 1:ti){

              Vtilde00 = matrix(Vtilde[Tinf0==ti,(1:ti)],sum(Tinf0==ti),ti)
              deriv_C_S[Tinf0==ti] = deriv_C_S[Tinf0==ti] + Xtilde[Tinf0==ti,t]*Vtilde[Tinf0==ti,t]*C_S_fun(S[Tinf0==ti]-1,
                                                                                                            matrix(Vtilde00[,(1:ti)!=t],sum(Tinf0==ti), sum((1:ti)!=t)  )  );
            }
          }
        }
      }else{
        # ti  =3
        # t=2
        for(ti in grid_T1){
          if(ti==1){
            deriv_C_S[Tinf0==ti] = deriv_C_S[Tinf0==ti] + Xtilde[Tinf0==ti,t,k]*Vtilde[Tinf0==ti,t];
          }else{
            for( t in 1:ti){

              Vtilde00 = matrix(Vtilde[Tinf0==ti,1:ti],sum(Tinf0==ti),ti)
              deriv_C_S[Tinf0==ti] = deriv_C_S[Tinf0==ti] + Xtilde[Tinf0==ti,t,k]*Vtilde[Tinf0==ti,t]*C_S_fun(S[Tinf0==ti]-1,
                                                                                                              matrix(Vtilde00[,(1:ti)!=t],sum(Tinf0==ti), sum((1:ti)!=t)  )  );
            }
          }
        }
      }

      # C_S_fun(t-1, Vtilde0[,(1:(T-1))!=u]-1)
      deriv_e = zeros(n,Tmax+1);
      if(length(dim( Xtilde)) == 2){
        for(ti in grid_T1){
          if(ti!=1){
            # Vtilde0 = as.matrix(Vtilde[Tinf0==ti,-c(ti)])
            Vtilde0 =matrix(Vtilde_min1 [Tinf0==ti,],sum(Tinf0==ti),dim(Vtilde_min1)[2])
            for( t in 1:(ti-1)){
              for( u in 1:(ti-1)){

                Vtilde00 = matrix(Vtilde0[,(1:(ti-1))],dim(Vtilde0)[1],ti-1)
                deriv_e[Tinf0==ti,t+1] = deriv_e[Tinf0==ti,t+1] + Xtilde[Tinf0==ti,u]*Vtilde[Tinf0==ti,u]* C_S_fun(t-1, matrix(Vtilde00[,(1:(ti-1))!=u],dim(Vtilde00)[1],sum((1:(ti-1))!=u)));
              }
            }
          }
        }

      }else{
        for(ti in grid_T1){
          if(ti!=1){
            # Vtilde0 = as.matrix(Vtilde[Tinf0==ti,-c(ti)])
            Vtilde0 = matrix(Vtilde_min1 [Tinf0==ti,],sum(Tinf0==ti),dim(Vtilde_min1)[2])
            for( t in 1:(ti-1)){
              for( u in 1:(ti-1)){

                Vtilde00 = matrix(Vtilde0[,(1:(ti-1))],dim(Vtilde0)[1],ti-1)
                deriv_e[Tinf0==ti,t+1] = deriv_e[Tinf0==ti,t+1] + Xtilde[Tinf0==ti,u,k]*Vtilde[Tinf0==ti,u]* C_S_fun(t-1, matrix(Vtilde00[,(1:(ti-1))!=u],dim(Vtilde00)[1],sum((1:(ti-1))!=u)));
              }
            }
          }
        }
      }
      deriv_lambda = zeros(n,Tmax+2);
      for(ti in grid_T1){
        deriv_lambda[Tinf0==ti,3:(ti+2)] = deriv_e[Tinf0==ti,2:(ti+1)]-deriv_e[Tinf0==ti,1:ti];
      }

      # dim(deriv_e)
      # Derivative of the coeff a_t wrt beta
      deriv_a = deriv_lambda[,1:(Tmax+1)]
      for(ti in grid_T1){
        deriv_a[Tinf0==ti,1:(ti+1)] = deriv_lambda[Tinf0==ti,1:(ti+1)] - matrix(deriv_lambda[Tinf0==ti,(ti+2)])%*%cheb_coeff[[ti]][1:(ti+1)]/(2^(ti+1));
      }
      # moy_deriv = (T+1)*mean(mat_combin*(deriv_a - bsxfun(@times, Mat_a, deriv_C_S/C_S_vec)),'all');
      moy_deriv = sum(mat_combin*(deriv_a -  Mat_a* ((deriv_C_S/C_S_vec)%*%rep(1,dim(Mat_a)[2]))),na.rm=T)/n;


      for_term2 = for_term2 + moy_deriv*  phi_b [!is.na(Tall),k];
      # term2 = b_hat * moy_deriv * phi;

    } ## end of the for k loop

    if(dimX==1){
      term2 = for_term2 * b_hat
    }else{
      term2 = for_term2 %*% b_hat
    }

    infl =  term1 + term2


    nb_clust =length(unique(Call1))
    et = matrix(0,1,dimX)
    for(ci in 1:nb_clust){
      if(dimX==1){
        et =et + var(infl[Call10==ci ,])
      }else{
        et =et + diag(var(infl[Call10==ci ,]))
      }
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

    # Delta_hat - length_CI/2
    # Save results


  }else{ ### end of Option 1 "quick".

    n <- dim(Yall)[1]
    nb_clust =length(unique(Call1))

    boundsall = vector("list")
    infl = vector("list")
    for(k in 1:dimX){
      infl[[k]]=matrix(0,n,2)
      # phi_b[[k]] = matrix(0,n,2)
    }

    # phi_b =matrix(0,n,dimX)

    sfInit(parallel = TRUE, cpus = nbCores, type = "SOCK")
    sfExportAll()
    sfLibrary(MarginalFElogit)
    sfLibrary(R.matlab)
    sfLibrary(pracma)

    Tall0 = Tall[!is.na(Tall)]
    grid_T1 = grid_T[grid_T >=min( Tall0 )]
    n_s <- sum(!is.na(Tall))

    # X = Xall[!is.na(Tall),,]
    # Y = matrix(Yall[!is.na(Tall),], sum(!is.na(Tall)), dim(Yall)[2])
    # Tinf0 = Tinf[!is.na(Tall)]
    # Call10 = matrix(Call1[!is.na(Tall),], sum(!is.na(Tall)), 1)


    for_average = vector("list")

    #
    # ti0=1
    # ci =1
    count = 1
    for(ci in 1:nb_clust){
      for(ti0 in 1:length(grid_T1)){

        T = grid_T1[ti0]

        selecti = Tinf==T & Call1==ci & !is.na(Tall)

        Y = matrix(Yall[selecti,1:T], sum(selecti) , T )
        X = array(Xall[selecti, 1:T ,], c(sum(selecti), T ,dimX))
        Tall0 <- Tall[selecti]

        # % Compute the influence function of beta_hat. Useful for inference on
        # % Delta, at the end.
        phi  = phi_b[selecti,]

        n <- dim(X)[1]
        T <- dim(X)[2]


        if(length(dim(X))==3){
          dimX = dim(X)[3]
        }else{
          dimX=1
        }

        T_minus_t = repmat(T - (0:T),n,1);
        S = apply(Y,1,sum)

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

        #
        # # % Parameter used to compute numerical derivatives wrt beta and gamma
        hderiv=1e-3;
        #

        # % Initialisation
        C_mat = zeros(n,T+1);
        PSt_X = zeros(n,T+1);
        # % Step 1: CMLE
        # options=optimoptions('fminunc','Display','off');
        # b_hat = optim(par = beta0, log_lik_FE_logit,Y=Y,X=X)$par


        index = zeros(n,T)

        # % Step 2: bounds on the AME

        # % We first need nonparametric estimators of mt(x).
        # % We use mt(x)=ct(x)/c0(x) for t\geq 0, with (c0(x),...,cT(x))'
        # % = A x (Pr(S=0|X=x)/C0(x;beta0),...,Pr(S=T|X=x)/CT(x;beta0).

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

        # %% Computation of the bandwidth for local linear estimation.
        # % We fix h so that estimated avg absolute bias = ratio * estimated
        # % avg asymp. std dev. To estimate both, we assume that the fixed
        # % effect is actually constant (and then estimated by intercept).


        # % Uniform draws on Supp(X). Used below to get the bandwidth.

        n_mx = 2000
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
        index_bis= matrix(0,n_mx,T)
        for(k in 1:dimX){
          index_bis = index_bis + matrix(randunif_x[,,k] * b_hat[k],n_mx,T);
        }

        # t = 2
        # ratio = 1/100
        h = zeros(1,T+1);
        for (t in 0:T){
          C_mat[,t+1] = C_S_fun(t,V);
          # % Computation of the avg std dev: (R(k)^q int sigma^2(x)dx/n)^{1/2}
          # % (see Hansen's notes). We use the uniform draws above. The
          # % estimator is an average of sigma^2(x) over these draws
          mean_approx_St = exp(t*intercept)*apply(1-Lambda(intercept+index_bis),1,prod)* C_S_fun(t,exp(index_bis));
          std_approx_St = scaling0*sqrt(RK^(dimX*T) * mean(mean_approx_St*(1-mean_approx_St))/n);

          # % Computation of the average absolute bias
          mean_approx_St = exp(t*intercept)*apply(Lambda(index_int),1,prod)*C_mat[,t+1];
          if (t==0){
            deriv_log_C = 0;
          }else{
            bbx = V* (matrix(C_mat[,t])%*%rep(1,dim(V)[2])) #bsxfun(@times, V,  C_mat[,t])
            deriv_log_C = b_hat * ( bbx / (matrix(C_mat[,t+1])%*%rep(1,dim(bbx)[2])))  # bsxfun(@rdivide, ,C_mat[,t+1]);
          }
          term1 = (deriv_log_C + b_hat * Lambda(index_int))^2;
          term2 = deriv_log_C*(b_hat - deriv_log_C) + b_hat^2*Lambda(index_int,1);
          bias_approx_St = sum(mean(abs( (as.matrix(mean_approx_St)%*%rep(1,dim(term1)[2]))*(term1+term2))));

          # % Choice of the bandwidth
          h[1,t+1] = (std_approx_St / (sqrt(ratio)*bias_approx_St))^(1/(2+dimX*T/2));
          # h[1,t+1] = (ratio*std_approx_St / bias_approx_St)^(1/(2+dimX*T/2));
        }


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

        ### case PSt_X==0
        temp0 = rowSums(abs(PSt_X))
        totreat = (1:dim(PSt_X)[1])[temp0==0]

        if(length(totreat)>0){

          # i0=1
          for(i0 in 1:length(totreat)){
            i1 = totreat[i0]
            values = matrix(0,dim(X)[1], 1)
            for(k in 1:dimX){
              values =  values + rowSums((X[,,k]-matrix(1,dim(X)[1],1)%*%X[i1,,k])^2)
            }

            nonzero <- cbind(values,!(temp0==0), 1:dim(PSt_X)[1])
            nonzero0 <- nonzero [nonzero[,2]==1,]
            ref = nonzero0 [which.min(nonzero0 [,1]),]
            i00 = ref[3]

            PSt_X[i1,] = PSt_X[i00,]
          }
        }

        # % Normalization: to be sure that the probabilities sum to 1.
        # PSt_X = bsxfun(@rdivide, PSt_X, sum(PSt_X,2));
        PSt_X = PSt_X / (matrix(apply(PSt_X,1,sum))%*%rep(1,dim(PSt_X)[2]));

        T_minus_t = repmat(T - (0:T),n,1);


        # deriv_gamma_all = vector("list")
        # version ="new"
        # s=1
        # if(version == "new"){


        # # % Computation of the bounds and Pr(Ihat(m_i)=I(m_i))
        out =  Delta_bounds(b_hat, V, S, dimX, PSt_X, f, h,mat_C_k_T[[T]], RK, T_minus_t,Tall0,X)
        boundsall[[count]] = out[[1]]
        bounds = out[[1]]
        bounds_ind = out[[2]]
        c_hat = out[[3]]


        # % Computation of the influence function
        Z_EZX = matrix(S)%*%rep(1,T+1)==(matrix(rep(1,length(S)))%*%t(matrix(0:T))) -  PSt_X;

        # k=1
        for(k in 1:dimX){
          deriv_beta = zeros(dimX,2);
          # i=1
          for (i in 1:dimX){
            out = Delta_bounds(b_hat+hderiv*((1:dimX)==i),V, S, dimX, PSt_X, f, h,mat_C_k_T[[T]], RK, T_minus_t,Tall0,X);
            for_deriv_b= out[[1]]
            # bounds_ind = out[[2]]
            # tx_class[s] = out[[3]]
            deriv_beta[i,]= (for_deriv_b[k,] - bounds[k,])/hderiv;
          }

          deriv_gamma = zeros(T+1,2)

          for (t in 0:T){
            out = Delta_bounds(b_hat, V, S, dimX, PSt_X+hderiv*ones(n,1)%*%((0:T)==t), f, h, mat_C_k_T[[T]], RK, T_minus_t,Tall0,X)
            for_deriv_g <- out[[1]]
            deriv_gamma[t+1,] = (for_deriv_g[k,] - bounds[k,] )/hderiv;
          }

          # % Term corresponding to the influence of the nonparametric
          # % estimation of the Pr(S=t|X)
          # if(is.null(infl[k])){
          # infl[[k]] = rbind( bounds_ind[,k] + phi %*% deriv_beta + Z_EZX %*% deriv_gamma)
          # }else{
          infl[[k]][selecti ,] =  bounds_ind[,k] + phi %*% deriv_beta + Z_EZX %*% deriv_gamma
          # dim(bounds_ind[,k] + phi %*% deriv_beta + Z_EZX %*% deriv_gamma)
          # dim( infl[[k]])
          # }
        }
        count=count+1




      } ### end of T_i loop
    }

    sfStop()

    count=1
    Delta_hat= matrix(0, nb_var,2)
    for(ci in 1:nb_clust){
      for(ti in grid_T1){
        Delta_hat=  Delta_hat+ boundsall[[count]][selectX,]*sum(selecti)/n_s
        count=count+1
      }
    }


    for_average = vector("list")
    count=1
    Delta_hat0= matrix(0, nb_var,2)
    for(ci in 1:nb_clust){
      for(ti in grid_T1){
        Delta_hat0=  Delta_hat0 + boundsall[[count]][selectX,]*sum(selecti)/(n_s*ti)
        count=count+1
      }
    }
    for_average[[1]] <- Delta_hat0



    # Delta_hat[,1] - length_CI/2
    length_CI= matrix(NA,1, nb_var)
    et=matrix(NA,dimX,2)
    for(k in selectX){
      et[k,] =apply(infl[[k]][!is.na(Tall),],2,std)/sqrt(n_s);
    }
    et=matrix(et[selectX,], nb_var,2)
    bias_sup = matrix(0,2,1)

    # stop =  Sys.time()
    # time =  start_time -stop;

    CI = matrix(NA, nb_var,2)

    ## compute CI at level alpha.
    # k=1
    for(k in 1: nb_var){
      quant = quantile_IM(alpha, Delta_hat[k,], et[k,]);

      # phi_alpha=1{sqrt(n)|beta hat_k|/std(phi_k)>q_{1-alpha/2}}
      phi_alpha = sqrt(n_s)*abs(b_hat[k])/std(infl[[k]])> qnorm(1-alpha/2)
      if(phi_alpha){
        CI[k,] = cbind(Delta_hat[k,1] - quant* et[k,1], Delta_hat[k,2] + quant*et[k,2]);
      }else{
        CI[k,] = cbind(min(0,Delta_hat[k,1] - quant* et[k,1]), max(0,Delta_hat[k,2] + quant*et[k,2]));
      }

      phi_alpha_m[k] = phi_alpha*1
      length_CI[k] = CI[k,2] - CI[k,1]

    }

    std_b = apply(phi_b,2,std)/sqrt(dim(phi_b)[1])
    # CI = cbind(bounds[,1] - quant* std_bounds[,1], bounds[,2] + quant*std_bounds[,2]);

  }

  end_time <- Sys.time()
  time = end_time - start_time


  out <- vector("list")
  out[[1]] <- Option
  out[[2]] <- n_s
  out[[3]] <- grid_T1
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
  out[[15]] <-  c_hat

  names(out) <- c("Option","n","Tmax","Time","Delta_hat","length_CI","et","bias_sup","CI","phi_alpha","b_hat","std_b","influence","for_average","c_hat")

  return(out)
}
