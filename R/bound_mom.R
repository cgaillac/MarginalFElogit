#' Bound on the T0+1-th moment given the vector of first moments m.
#'
#' @param m the vector m until m_{T0}
#' @param T0 the lenght of the vector of moments m
#' @param indic We may have length(m)<T0 if m is at the boundary. In this case, this 3rd
#' argument (indic) is required, with indic=1 if the next element of m is =
#'    to the upper bound, indic=-1 if it is = to the lower bound
#'
#'
#' @return return the bounds on the T0+1-th moment
#' @export
#'
# @examples
bound_mom <- function(m,T0,indic){
  # m = matrix(mi,length(mi),1)
  # T0 = T

  res=matrix(0,2,1)
  out <- tryCatch({
    if (length(m)<T0){
      # % Boundary case: qinf=qsup
      # % We complete the vector m until m_{T+1}
      T1 = length(m);
      if (indic==-1){
        # % We are at the lower bound
        if (mod(T1,2)==0){
          # % Even case
          for (t in T1:T0){
            m0=m[(t-T1+1):(t-T1/2+1)]
            m1=m[(t-T1/2+1):length(m)];
            m10=c(m1,0);
            y1=det(hankel(m0,m10))
            m11=c(m1,1);
            y2=det(hankel(m0,m11));
            m[t+1] = y1/(y1-y2);
          }
        }else{
          # % Odd case
          for (t in T1:T0){
            if( t==T1){
              m0=c(1,m[1:((T1+1)/2)])
            }else{
              m0 = m[(t-T1):(t-(T1-1)/2)]
            }
            m1=m[(t-(T1-1)/2):length(m)];
            m10=c(m1,0);
            y1=det(hankel(m0,m10));
            m11=c(m1,1);
            y2=det(hankel(m0,m11));
            m[t+1] = y1/(y1-y2);
          }
        }
      }else{
        # % We are at the upper bound
        if (mod(T1,2)==0){
          # % Even case
          for (t in T1:T0){
            if (t==T1){
              m0=c(1,m[1:T1/2]) - m[1:(T1/2+1)]
            }else{
              m0=m[(t-T1):(t-T1/2)] - m[(t-T1+1):(t-T1/2+1)];
            }
            m1=m[(t-T1/2):(length(m)-1)] - m[(t-T1/2+1):length(m)];
            m10 = c(m1,m[length(m)]);
            y1=det(hankel(m0,m10));
            m11=c(m1,m[length(m)]-1);
            y2=det(hankel(m0,m11));
            m[t+1] = y1/(y1-y2);
          }
        }else{
          # % Odd case
          if(T1==1){
            #% Bernoulli case
            m=matrix(m[1],T+1,1)
          }else{
            for (t in T1:T0){
              m0=m[(t-T1+1):(t-(T1-1)/2)] - m[(t-T1+2):(t-(T1-3)/2)];
              m1=m[(t-(T1-1)/2):(length(m)-1)] - m[(t-(T1-3)/2):length(m)];
              m10 =c(m1,m[length(m)]);
              y1=det(hankel(m0,m10));
              m11=c(m1,m[length(m)]-1);
              y2=det(hankel(m0,m11));
              m[t+1] = y1/(y1-y2);
            }
          }
        }
      }
      res[1]=m[length(m)];
      res[2]=m[length(m)];
    }else{
      # % Interior case: qinf < qsup
      if(length(m)>T0){
        m=m[1:T0]
      }

      # % T even
      if (mod(T0,2)==0){
        # % Lower bound
        m0=m[1:(T0/2+1)]
        m1=m[(T0/2+1):length(m)]
        m10=c(m1,0)
        y1=det(hankel(m0,m10))
        m11=c(m1,1)
        y2=det(hankel(m0,m11))
        res[1] = y1/(y1-y2)

        m0=c(1,m[1:(T0/2)]) - m0;
        m1=m[(T0/2):(length(m)-1)] - m[(T0/2+1):length(m)];
        m10 = c(m1,m[length(m)]);
        y1=det(hankel(m0,m10));
        m11=c(m1,m[length(m)]-1);
        y2=det(hankel(m0,m11));
        res[2] = y1/(y1-y2);
      }else{
        # % T odd
        m0=c(1,m[1:((T0+1)/2)])
        m1=m[((T0+1)/2):length(m)]
        m10=c(m1,0)
        y1=det(hankel(m0,m10))
        m11=c(m1,1)
        y2=det(hankel(m0,m11));
        res[1]= y1/(y1-y2)

        if(T0==1){
          y1=m[1]
          y2=m[1] -1;
          res[2] = y1/(y1-y2);

        }else{
          m0=m[1:((T0+1)/2)]- m[2:((T0+3)/2)]
          m1=m[((T0+1)/2):(length(m)-1)] - m[((T0+3)/2):length(m)];
          m10 =c(m1,m[length(m)]);
          y1=det(hankel(m0,m10));
          m11=c(m1,m[length(m)]-1);
          y2=det(hankel(m0,m11));
          res[2] = y1/(y1-y2);
        }
      }
    }
  },error=function(cond){

  })


  if(sum(is.na(res))>0){
    if(m[1]==0){
      res=matrix(0,2,1)
    }else{
      res=bound_mom(m[1],T,1);
    }
  }
  return(res)
}
