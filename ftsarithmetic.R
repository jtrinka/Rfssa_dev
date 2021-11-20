#' Addition of Functional Time Series
#'
#' A method that lets you perform functional time series (\code{\link{fts}}) addition and scalar addition.
#' @return an object of class \code{\link{fts}}.
#'
#' @param Y1 an object of class \code{\link{fts}} or scalar
#' @param Y2 an object of class \code{\link{fts}} or scalar
#' @examples
#'
#' \dontrun{
#' require(fda)
#' require(Rfssa)
#' data(Callcenter) # Read data
#' u=seq(0,1,length.out=240) # Define domain of functional data
#' d=12 # number of basis elements
#' basis=create.bspline.basis(rangeval = c(0,1),nbasis = d) # create basis object
#' smooth.calls=smooth.basis(u, matrix(nrow=240,ncol=365,Callcenter$calls), basis)
#' Y=fts(smooth.calls$fd) # create functional time series
#' plot(Y)
#' Yplus=Y+Y # add the functional time series to itself
#' plot(Yplus)
#' Yplus2=Y+2 # add 2 to every term in the functional time series
#' plot(Yplus2)
#' }
#'
#' @seealso \code{\link{fts}}
#' @useDynLib Rfssa
#' @export
'+.fts' <- function(Y1,Y2){


  if(class(Y1)=="fts" && class(Y2)=="fts"){

    if(length(Y1@C)!=length(Y2@C)){

      stop("Functional time series are of different length.")

    }


    p=length(Y1@C)
    grid=Y1@grid
    basis=Y1@B
    eval_fts=list()
    for(j in 1:p){
      if(max(max(Y1@B[[j]]-Y2@B[[j]]))>=10^-14){
        stop(paste("The basis matrix corresponding to variable ",as.character(j)," should be the same between fts objects.",sep=""))
      }

      eval_fts[[j]]=basis[[j]]%*%(Y1@C[[j]]+Y2@C[[j]])


    }
  }else if(class(Y1)=="fts" && class(Y2)=="numeric" || class(Y1)=="fts" && class(Y2)=="matrix"){

    Y2=as.numeric(Y2)
    p=length(Y1@C)
    if(length(Y2)==1) Y2=rep(Y2,p);
    grid=Y1@grid
    basis=Y1@B
    eval_fts=list()
    for(j in 1:p){

      eval_fts[[j]]=basis[[j]]%*%(Y1@C[[j]])+Y2[j]

    }


    }else if(class(Y1)=="numeric" && class(Y2)=="fts" || class(Y1)=="matrix" && class(Y2)=="fts"){

      Y1=as.numeric(Y1)
      p=length(Y2@C)
      if(length(Y1)==1) Y1=rep(Y1,p);
      grid=Y2@grid
      basis=Y2@B
      eval_fts=list()
      for(j in 1:p){

        eval_fts[[j]]=basis[[j]]%*%(Y2@C[[j]])+Y1[j]


    }


  }

  out=Rfssa::fts(eval_fts,basis,grid)
  return(out)

}

#' Subtraction of Functional Time Series
#'
#' A method that lets you perform functional time series (\code{\link{fts}}) subtraction and scalar subtraction.
#' @return an object of class \code{\link{fts}}.
#'
#'
#' @param Y1 an object of class \code{\link{fts}} or scalar
#' @param Y2 an object of class \code{\link{fts}} or scalar
#' @examples
#'
#' \dontrun{
#' require(fda)
#' require(Rfssa)
#' data(Callcenter) # Read data
#' u=seq(0,1,length.out=240) # Define domain of functional data
#' d=12 # number of basis elements
#' basis=create.bspline.basis(rangeval = c(0,1),nbasis = d) # create basis object
#' smooth.calls=smooth.basis(u, matrix(nrow=240,ncol=365,Callcenter$calls), basis)
#' Y=fts(smooth.calls$fd) # create functional time series
#' plot(Y)
#' Yminus=Y[4:8]-Y[14:18] # subtract elements of the functional time series from each other
#' plot(Yminus)
#' Yminus2=Y-2 # add 2 to every term in the functional time series
#' plot(Yminus2)
#' }
#'
#' @seealso \code{\link{fts}}
#' @useDynLib Rfssa
#' @export
'-.fts' <- function(Y1,Y2){


  if(class(Y1)=="fts" && class(Y2)=="fts"){

    if(length(Y1@C)!=length(Y2@C)){

      stop("Functional time series are of different length.")

    }

    p=length(Y1@C)
    grid=Y1@grid
    basis=Y1@B
    eval_fts=list()
    for(j in 1:p){
      if(max(max(Y1@B[[j]]-Y2@B[[j]]))>=10^-14){
        stop(paste("The basis matrix corresponding to variable ",as.character(j)," should be the same between fts objects.",sep=""))
      }

      eval_fts[[j]]=basis[[j]]%*%(Y1@C[[j]]-Y2@C[[j]])


    }
  }else if(class(Y1)=="fts" && class(Y2)=="numeric" || class(Y1)=="fts" && class(Y2)=="matrix"){

    Y2=as.numeric(Y2)
    p=length(Y1@C)
    if(length(Y2)==1) Y2=rep(Y2,p);
    grid=Y1@grid
    basis=Y1@B
    eval_fts=list()
    for(j in 1:p){

      eval_fts[[j]]=basis[[j]]%*%(Y1@C[[j]])-Y2[j]

    }


  }else if(class(Y1)=="numeric" && class(Y2)=="fts" || class(Y1)=="matrix" && class(Y2)=="fts"){

    Y1=as.numeric(Y1)
    p=length(Y2@C)
    if(length(Y1)==1) Y2=rep(Y1,p);
    grid=Y2@grid
    basis=Y2@B
    eval_fts=list()
    for(j in 1:p){

      eval_fts[[j]]=Y1[j]-basis[[j]]%*%(Y2@C[[j]])


    }


  }

  out=Rfssa::fts(eval_fts,basis,grid)
  return(out)

}


#' Multiplication of Functional Time Series
#'
#' A method that lets you multiply functional time series (\code{\link{fts}}) and perform scalar multiplication of functional time series.
#' @return an object of class \code{\link{fts}}
#'
#'
#' @param Y1 an object of class \code{\link{fts}} or scalar
#' @param Y2 an object of class \code{\link{fts}} or scalar
#' @examples
#'
#' \dontrun{
#' require(fda)
#' require(Rfssa)
#' data(Callcenter) # Read data
#' u=seq(0,1,length.out=240) # Define domain of functional data
#' d=12 # number of basis elements
#' basis=create.bspline.basis(rangeval = c(0,1),nbasis = d) # create basis object
#' smooth.calls=smooth.basis(u, matrix(nrow=240,ncol=365,Callcenter$calls), basis)
#' Y=fts(smooth.calls$fd) # create functional time series
#' plot(Y)
#' Ytimes=Y*Y # elementwise multiplication of the functional time series with itself
#' plot(Ytimes)
#' Ytimes2=2*Y # multiply 2 with every term in the functional time series
#' plot(Ytimes2)
#' }
#'
#' @seealso \code{\link{fts}}
#' @useDynLib Rfssa
#' @export
'*.fts' <- function(Y1,Y2){


  if(class(Y1)=="fts" && class(Y2)=="fts"){

    if(length(Y1@C)!=length(Y2@C)){

      stop("Functional time series are of different length.")

    }

    p=length(Y1@C)
    grid=Y1@grid
    basis=Y1@B
    eval_fts=list()
    for(j in 1:p){

      if(max(max(Y1@B[[j]]-Y2@B[[j]]))>=10^-14){
        stop(paste("The basis matrix corresponding to variable ",as.character(j)," should be the same between fts objects.",sep=""))
      }

      d = ncol(basis[[j]])
      Y1_j_eval=basis[[j]]%*%Y1@C[[j]]
      Y2_j_eval=basis[[j]]%*%Y2@C[[j]]
      eval_fts[[j]]=Y1_j_eval*Y2_j_eval


    }

  }else if(class(Y1)=="fts" && class(Y2)=="numeric" || class(Y1)=="fts" && class(Y2)=="matrix"){

    Y2=as.numeric(Y2)
    p=length(Y1@C)
    if(length(Y2)==1) Y2=rep(Y2,p);
    grid=Y1@grid
    basis=Y1@B
    eval_fts=list()
    for(j in 1:p){
      d = ncol(basis[[j]])
      Y1_j_eval=basis[[j]]%*%Y1@C[[j]]
      eval_fts[[j]]=Y1_j_eval*Y2[j]


    }


  }else if(class(Y1)=="numeric" && class(Y2)=="fts" || class(Y1)=="matrix" && class(Y2)=="fts"){

    Y1=as.numeric(Y1)
    p=length(Y2@C)
    if(length(Y1)==1) Y1=rep(Y1,p);
    grid=Y2@grid
    basis=Y2@B
    eval_fts=list()
    for(j in 1:p){
      d = ncol(basis[[j]])
      Y2_j_eval=basis[[j]]%*%Y2@C[[j]]
      eval_fts[[j]]=Y2_j_eval*Y1[j]


    }


  }

    out=Rfssa::fts(eval_fts,basis,grid)
    return(out)

}


#' Indexing into Functional Time Series
#'
#' A method that lets you index into a functional time series (\code{\link{fts}}).
#' @return an object of class \code{\link{fts}}
#'
#'
#' @param Y an object of class \code{\link{fts}}
#' @param i index
#' @examples
#'
#' \dontrun{
#' require(fda)
#' require(Rfssa)
#' data(Callcenter) # Read data
#' u=seq(0,1,length.out=240) # Define domain of functional data
#' d=12 # number of basis elements
#' basis=create.bspline.basis(rangeval = c(0,1),nbasis = d) # create basis object
#' smooth.calls=smooth.basis(u, matrix(nrow=240,ncol=365,Callcenter$calls), basis)
#' Y=fts(smooth.calls$fd) # create functional time series
#' Yind=Y[4:8] # take only the 4th through 8th functions
#' plot(Yind)
#' Yminus=Y[4:8]-Y[14:18] # subtract functions from each other
#' plot(Yminus)
#' }
#'
#'
#'
#' @seealso \code{\link{fts}}
#' @useDynLib Rfssa
#' @note can use ':' as an operator to specify a range of indices
#' @export
'[.fts'<-function(Y,i='index'){
  p=length(Y@C)
  grid=Y@grid
  basis=Y@B
  eval_fts=list()
  for(j in 1:p){

    eval_fts[[j]]=basis[[j]]%*%Y@C[[j]][,i]


  }

  out=Rfssa::fts(eval_fts,basis,grid)
  return(out)

}


