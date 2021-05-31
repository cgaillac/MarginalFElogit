#' This function truncates m_hat until it belongs to the moment space. The
#' truncated moment is put in m_trun.
#' Sigma is an estimator of the "variance" of m_hat, or more accurately, of
#' the asymptotic variance of m_hat, properly normalized.
#'
#' @param m_hat the estimated sequence of moments
#' @param Sigma is an estimator of the "variance" of m_hat, or more accurately, of the asymptotic variance of m_hat, properly normalized.
#' @param n the sample size
#' @param ind_Sig_pos indicates if Sigma is invertible
#'
#'  @return A list containing:
#'
#'    - m_trun: the truncated value of m_hat, which belongs to the moment space;
#'
#'    - indic=0 if m_hat was not truncated. Otherwise, indic=1 if we are above the upper bound, indic=-1 if we are  below the lower bound.
#'
#' @export
#'
# @examples
trunc_m <- function (m_hat, Sigma =NULL, n=NULL,ind_Sig_pos = NULL){
  # m_hat = matrix(mi)
  # Sigma = Sigma_m
  out = NULL
  out1 <- tryCatch({

    T=length(m_hat);

    if( !is.null(n)){
      ct = sqrt(2*log(log(n)));
    }

    k=0
    indic = 0
    m0 = rbind(1,m_hat)

    if( is.null(ind_Sig_pos)){
      if( is.null(Sigma) && is.null(n)){
        Sigma = 0;
        ind_Sig_pos = FALSE
      }else{
        ind_Sig_pos = det(Sigma)>0
      }
    }else{
      ind_Sig_pos = FALSE
    }

    while( (indic == 0) && (k < T)){
      k=k+1;
      ell = round(k/2);
      if( rem(k,2)==1){
        Bkp1 = hankel(m_hat[1:ell],m_hat[ell:k]);
        Akp1 = hankel(m0[1:ell],m0[ell:k]);
        Hinf = det(Bkp1);
        Hsup = det(Akp1 - Bkp1);

        if( !is.null(Sigma) && !is.null(n) && ind_Sig_pos){
          # % We apply the Jacobi formula to estimate the variance of Hinf and
          # % Hsup
          L = dim(Bkp1)[1];
          Jac_inf = zeros(1,T);
          Adj_inf = adjoint1(Bkp1);
          # i =1
          for( i in seq(1,T)){
            bx = matrix(1:L)%*%rep(1,L) + matrix(rep(1,L))%*%(1:L)
            Jac_inf[i] = sum(diag(Adj_inf * (bx ==(i+1))));
          }
          s_inf = sqrt(Jac_inf%*% Sigma%*% t(Jac_inf));
          Jac_sup = zeros(1,T);
          Adj_sup = adjoint1(Akp1 - Bkp1);
          for (i in seq(1,T)){
               bx = matrix(1:L)%*%rep(1,L) + matrix(rep(1,L))%*%(1:L)
              Jac_sup[i] = sum(diag(Adj_sup * ((bx ==(i+2)) - (bx ==(i+1) ))));
          }
          s_sup = sqrt(Jac_sup%*% Sigma %*% t(Jac_sup));
          if( min(Hinf/s_inf,Hsup/s_sup)<=ct){
                  indic = 2*(Hinf/s_inf>Hsup/s_sup) - 1;
           }
          }else{
              if(min(Hinf,Hsup)<=1e-10){
                  indic = 2*(Hinf>Hsup) - 1}
          }
        }else{
          Akp1 = hankel(m0[1:(ell+1)],m0[(ell+1):(k+1)]);
          Bk = hankel(m_hat[1:ell],m_hat[ell:(k-1)]);
          Ck = hankel(m_hat[2:(ell+1)],m_hat[(ell+1):k])
          Hinf = det(Akp1);
          Hsup = det(Bk - Ck);

        if(  !is.null(Sigma) && !is.null(n) && ind_Sig_pos){
              L = dim(Akp1)[1];
              Jac_inf = zeros(1,T);
              Adj_inf = adjoint1(Akp1);
              for( i in seq(1,T)){
                   bx = matrix(1:L)%*%rep(1,L) + matrix(rep(1,L))%*%(1:L)
                  Jac_inf[i] = sum(diag((Adj_inf * (bx ==(i+2)))));
              }
              s_inf = sqrt(Jac_inf %*% Sigma %*% t(Jac_inf));
              Jac_sup = zeros(1,T);
              Adj_sup = adjoint1(Bk - Ck);
              for( i in seq(1,T)){
                  # Jac_sup(i) = trace(Adj_sup * ((bsxfun(@plus,(1:L-1),...
                  #              (1:L-1)')==i+1) - (bsxfun(@plus,(1:L-1),...
                  #                                        (1:L-1)')==i)));
                bx = matrix(1:(L-1))%*%rep(1,(L-1)) + matrix(rep(1,L-1))%*%(1:(L-1))
                Jac_sup[i] = sum(diag(Adj_sup * ((bx ==(i+1)) - (bx ==(i) ))));
               }
              s_sup = sqrt(Jac_sup %*% Sigma %*% t(Jac_sup));
              if(min(Hinf/s_inf,Hsup/s_sup)<=ct){
                  indic = 2*(Hinf/s_inf>Hsup/s_sup) - 1;
              }

        }else if( min(Hinf,Hsup)<=1e-10){
                  indic = 2*(Hinf>Hsup) - 1;
        }
      }
    }

    if( k==1 & T>1){
      # % In this very peculiar case, the distribution is either the Dirac at 0
      #% or 1
      m_trun = (norm(m_hat-1, "I") < norm(m_hat, "I"))*1;
    }else if(indic!=0){
      m_trun = m_hat[1:(k-1)];
    }else{
      m_trun = m_hat;
    }

    out= vector("list")
    out[[1]] <- m_trun
    out[[2]] <- indic
  },
  error=function(cond) {

  })

  if(is.null(out)){
    out <- trunc_m(m_hat, Sigma =NULL, n=NULL, ind_Sig_pos = FALSE)
  }

return(out)
}
