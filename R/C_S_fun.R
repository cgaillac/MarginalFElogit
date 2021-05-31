#'Computes by induction one of the elementary symm. polynomial on each row
#'of V1, matrix of size n x T.
#'
#'If S1 is scalar, computes the S1-th elem. symm. polynomial
#'If S1 is a n-column vector, computes, for row i, the S1(i)-th elem. symm.
#'polynomial.
#'
#' @param S1 the vector of values for the statistic S.
#' @param V1 the matrix of size n x T which coefficients are
#'           used in the elementary symm. polynomial.
#'
#' @return a vector if S1 is scalar, or a matrix otherwise, containing
#' the elementary symm. polynomial on each row of V1.
#'
#' @export
#'
# @examples
C_S_fun <- function(S1, V1){

n=dim(V1)[1]
T = dim(V1)[2]

if (length(S1)==1){
  if (S1<0 || S1>T){
    res=0;
  }else if(S1==0){
    res=1;
  }else if(S1==T){
    if(dim(V1)[2]>1){
      res=matrix(apply(V1,1,prod),dim(V1)[1],1)
    }else{
      res=V1;
    }
  }else{
    # if(dim(V1)[2]<=2 || n==1){
    #    if(dim(V1)[2] ==1){
    #      VT_1 = matrix(0,1)
    #    }else{
    #      VT_1 = matrix(V1[,1:(dim(V1)[2]-1)],dim(V1)[1],dim(V1)[2]-1)
    #    }
    #     res=V1[,T]*C_S_fun(S1-1,VT_1)+C_S_fun(S1,VT_1);
    # }else{
    #   VT_1 = V1[,1:(dim(V1)[2]-1)]
    #   res=V1[,T]*C_S_fun(S1-1,VT_1)+C_S_fun(S1,VT_1);
    # }

    VT_1 = matrix(V1[,1:(dim(V1)[2]-1)],dim(V1)[1],dim(V1)[2]-1)
    #VT_1 = V1[,1:(dim(V1)[2]-1)]
    res=V1[,T]*C_S_fun(S1-1,VT_1)+C_S_fun(S1,VT_1);

  }
}else{
  res=zeros(n,1);
  if(sum(S1==0)!=0){
    res[S1==0]=1;
  }
  if(sum(S1==T)!=0){
    if(dim(V1)[2]>1){
      res[S1==T]=apply(matrix(V1[S1==T,],sum(S1==T),dim(V1)[2])   ,1,prod);
    }else{
      res[S1==T]=V1[(S1==T),];
    }
  }
  indic = (S1>0 & S1<T)
  if(sum(indic)>0){
    VT = matrix(V1[indic,T], sum(indic),1);
    # if(dim(V1)[2]<=2  || n==1){
    #   VT_1= matrix(V1[indic,1:(dim(V1)[2]-1)], sum(indic),dim(V1)[2]-1);
    # }else{
    #   VT_1= matrix( V1[indic,1:(dim(V1)[2]-1)], sum(indic), (dim(V1)[2]-1) )
    # }
    VT_1= matrix( V1[indic,1:(dim(V1)[2]-1)], sum(indic), (dim(V1)[2]-1) )
    Sind= matrix(S1[indic], sum(indic),1)
    res[indic]=VT*C_S_fun(Sind-1,VT_1)+C_S_fun(Sind,VT_1);
  }
}
return(res)
}
