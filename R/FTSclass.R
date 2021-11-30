#' Functional Time Series Class
#'
#' This function is used to create functional time series (fts) objects from lists of discretely sampled data, basis specifications, and grid elements which provide the domain that each variable is observed over. Each variable is assumed to be observed over a regular and equidistant grid.
#' @param X  A list of length p where p is a positive integer specifying the number of variables observed in the fts. Each entry in the list should be a matrix or
#' an array. For a variable observed over a one-dimensional domain, the list entry should be a m by N matrix where m is the number of sampling points for each observation and N is the number of observations.
#' For a variable observed over a two-dimensional domain, the list entry should be a m_1 by m_2 by N array where m_1 is the number of sampling points in the horizontal direction and
#' m_2 is the number of sampling points in the vertical direction.
#' @param B A list of length p. Each entry in the list should be either a matrix specifying the basis for each variable or each list entry should be a list specifying the number of basis elements and desired basis type to be used in the smoothing process.
#' For a variable observed over a one-dimensional domain, the list entry should be a m by d matrix or the list entry should be a list of the form, \code{list(d, "basistype")}, where d specifies the number basis elements and "basistype" is a string that is set to \code{"bspline"} or \code{"fourier"} to specify the type of basis.
#' For a variable observed over a two-dimensional domain, the list entry should be a (m_1)(m_2) by (d_1)(d_2) matrix or the list entry should be a list of the form, \code{list(d_1,d_2,"basistype1","basistype2")}, where d_1 specifies the number of basis elements in the horizontal direction, d_2 specifies the number of basis elements in the vertical direction, \code{"basistype1"} is set to be \code{"bspline"} or \code{"fourier"} which specifies the type of basis in the horizontal direction, and \code{"basistype2"} is set to be \code{"bspline"} or \code{"fourier"} which specifies the type of basis in the vertical direction.
#' @param grid A list of length p. Each entry in the list should either be a numeric or a list of numeric elements depending on the dimension of the domain the variable is observed over.
#' For a variable observed over a one-dimensional domain, the list entry should be either an ordered (from smallest to largest) numeric of length m giving the sampling points or a numeric of the form, \code{c(x_1,x_2)} where x_1 and x_2 are the smallest and largest values attained in the domain of the variable respectively. In addition, these list entries can also be provided in a list.
#' For a variable observed over a two-dimensional domain, the list entry should be either a list of the form, \code{list(u,v)}, or of the form, \code{list(c(x_1,x_2),v)}, or of the form, \code{list(u,c(x_3,x_4))}, or of the form, \code{list(c(x_1,x_2),c(x_3,x_4))} where u is a numeric with minimium and maximum values of x_1 and x_2 respectively and v is a numeric with minimium and maximum values of x_3 and x_4 respectively.
#' @importFrom fda fd inprod eval.fd smooth.basis is.fd create.bspline.basis create.fourier.basis eval.basis
#' @importFrom methods is new
#' @seealso \code{\link{fssa}}
#' @note Refer to \code{\link{fssa}} for an example on how to run this function.
#' @export
fts <- function(X, B, grid){
  fts=new(Class = "fts",C = X, B = B, grid = grid)
  p = length(X)
  # Loop over each variable in the fts
  for(j in 1L:p){
    # Deciding whether or not the variable corresponds to an fts observed over a two-dimensional domain or a one-dimensional domain.
    if(class(grid[[j]])=="list" && length(grid[[j]])==2){
      # Define a function used to build the grid for variables observed over a two-dimensional domain.

      twoDgrid=function(u,v){
        grid=matrix(nrow=(length(u)*length(v)),ncol=2,data=NA)
        for(i in 1:length(u)){
          for(j in 1:length(v)){
            k=(i-1)*length(v)+j
            grid[k,]=c(u[i],v[j])


          }
        }

        return(grid)

      }

      if(length(grid[[j]][[1]])==2 && length(grid[[j]][[2]])==2){
        # If the grid is specified with min and max values in both horizontal and vertical axes.
        # Example: list(c(1,2),c(3,4))

        u=seq(grid[[j]][[1]][1],grid[[j]][[1]][2],length.out=ncol(X[[j]]))
        v=seq(grid[[j]][[2]][1],grid[[j]][[2]][2],length.out=nrow(X[[j]]))
        M_x = length(u)
        M_y = length(v)
        x_min = min(u)
        x_max = max(u)
        y_min = min(v)
        y_max = max(v)

      }else if(length(grid[[j]][[1]])==2 && length(grid[[j]][[2]])>2){
        # If the grid is specified with min and max values in the horizontal axis and with the actual grid in the vertical axis.
        # Example: list(c(1,2),v)

        u=seq(grid[[j]][[1]][1],grid[[j]][[1]][2],length.out=ncol(X[[j]]))
        M_x = length(u)
        M_y = length(grid[[j]][[2]])
        x_min = min(u)
        x_max = max(u)
        y_min = min(grid[[j]][[2]])
        y_max = max(grid[[j]][[2]])
        v=grid[[j]][[2]]


      }else if(length(grid[[j]][[1]])>2 && length(grid[[j]][[2]])==2){
        # This case is similar as seen in the previous else if statement
        # Example: list(u,c(1,2))

        v=seq(grid[[j]][[2]][1],grid[[j]][[2]][2],length.out=nrow(X[[j]]))
        M_x = length(grid[[j]][[1]])
        M_y = length(v)
        x_min = min(grid[[j]][[1]])
        x_max = max(grid[[j]][[1]])
        y_min = min(v)
        y_max = max(v)
        u=grid[[j]][[1]]

      }else{

        # If the grid is specified with the actual grid in the horizontal and vertical axes.
        # Example: list(u,v)
        M_x = length(grid[[j]][[1]])
        M_y = length(grid[[j]][[2]])
        x_min = min(grid[[j]][[1]])
        x_max = max(grid[[j]][[1]])
        y_min = min(grid[[j]][[2]])
        y_max = max(grid[[j]][[2]])
        u=seq(x_min,x_max,length.out=M_x)
        v=seq(y_min,y_max,length.out=M_y)

      }
      # Create the grid using u and v
      fts@grid[[j]]=twoDgrid(u,v)

      # If the fts only contains one observation
      if(length(dim(X[[j]]))==2){

        X_temp=matrix(data=NA,nrow=M_x*M_y,ncol=1)
        count=1
        for(i in 1:M_x){
          for(k in 1:M_y){

            X_temp[count,1]=X[[j]][i,k]
            count=count+1

          }

        }
        X[[j]]=X_temp

      # If the fts contains more than one observation
      }else if(length(dim(X[[j]]))==3){

        N=dim(X[[j]])[3]
        X_temp=matrix(data=NA,nrow=M_x*M_y,ncol=N)
        for(n in 1:N){
          count=1
          for(i in 1:M_x){
            for(k in 1:M_y){

              X_temp[count,n]=X[[j]][i,k,n]
              count=count+1

            }

          }

        }

        X[[j]]=X_temp

      }
    # If the variable is observed over a one-dimensional domain
    }else if(class(grid[[j]])=="list" && length(grid[[j]])==1 && length(grid[[j]][[1]])>2){
      # If the grid is supplied in a list, the variable corresponds with a one-dimensional domain, and the list element is the grid.
      # Example: list(u)
      M_x=length(grid[[j]][[1]])
      u = seq(min(grid[[j]][[1]]),max(grid[[j]][[1]]),length.out=M_x)
      fts@grid[[j]]=matrix(u,ncol=1)
      X[[j]]=as.matrix(X[[j]])



    }else if(class(grid[[j]])=="list" && length(grid[[j]])==1 && length(grid[[j]][[1]])==2){

      # If the grid is supplied in a list, the variable corresponds with a one-dimensional domain, and the list element is the min and max values of the grid.
      # Example: list(c(1,2))
      u=seq(grid[[j]][[1]][1],grid[[j]][[1]][2],length.out=nrow(X[[j]]))
      fts@grid[[j]]=matrix(u,ncol=1)
      X[[j]]=as.matrix(X[[j]])



    }else if(class(grid[[j]])=="numeric" && length(grid[[j]])==2){

      # If the grid specified with the min and max values of the grid.
      # Example: c(1,2)
        u=seq(grid[[j]][1],grid[[j]][2],length.out=nrow(X[[j]]))
        fts@grid[[j]]=matrix(u,ncol=1)
        X[[j]]=as.matrix(X[[j]])


    }else{
      # If the grid specified with the grid itself.
      # Example: u
      M_x=length(grid[[1]])
      fts@grid[[j]]=matrix(grid[[j]],ncol=1)
      X[[j]]=as.matrix(X[[j]])


      }



    # Generating basis matrices
    if(class(B[[j]])=="list" && B[[j]][[2]]=="bspline" &&  ncol(fts@grid[[j]])==1){
      # Generating a bspline basis for variables whose domain is one-dimensional using a supplied list.
      # Example: list(10,"bspline")
      fts@B[[j]]=fda::eval.basis(evalarg = seq(min(fts@grid[[j]]),max(fts@grid[[j]]),length.out=nrow(fts@grid[[j]])),basisobj = fda::create.bspline.basis(rangeval = c(min(fts@grid[[j]]),max(fts@grid[[j]])),nbasis = B[[j]][[1]]))
      B[[j]]=fts@B[[j]]

    }else if(class(B[[j]])=="list" && B[[j]][[2]]=="fourier" &&  ncol(fts@grid[[j]])==1){
      # Generating a fourier basis for variables whose domain is one-dimensional using a supplied list.
      # Example: list(10,"fourier")
      fts@B[[j]]=fda::eval.basis(evalarg = seq(min(fts@grid[[j]]),max(fts@grid[[j]]),length.out=nrow(fts@grid[[j]])),basisobj = fda::create.fourier.basis(rangeval = c(min(fts@grid[[j]]),max(fts@grid[[j]])),nbasis = B[[j]][[1]]))
      B[[j]]=fts@B[[j]]
    }

    if(class(B[[j]])=="list" && B[[j]][[3]]=="bspline" && B[[j]][[4]]=="bspline" && length(grid[[j]])==2){
      # Generating a bspline basis for variables whose domain is two-dimensional using a supplied list.
      # Example: list(list(10,"bspline"),list(12,"bspline"))
      b_1=fda::eval.basis(evalarg = seq(min(fts@grid[[j]][,1]),max(fts@grid[[j]][,1]),length.out=M_x),basisobj = fda::create.bspline.basis(rangeval = c(min(fts@grid[[j]][,1]),max(fts@grid[[j]][,1])),nbasis = B[[j]][[1]]))
      b_2=fda::eval.basis(evalarg = seq(min(fts@grid[[j]][,2]),max(fts@grid[[j]][,2]),length.out=M_y),basisobj = fda::create.bspline.basis(rangeval = c(min(fts@grid[[j]][,2]),max(fts@grid[[j]][,2])),nbasis = B[[j]][[2]]))
      fts@B[[j]]=kronecker(b_1,b_2)
      B[[j]]=fts@B[[j]]

    }else if(class(B[[j]])=="list" && B[[j]][[3]]=="bspline" && B[[j]][[4]]=="fourier" && length(grid[[j]])==2){
      # Generating a bspline/fourier basis for variables whose domain is two-dimensional using a supplied list.
      # Example: list(list(10,"bspline"),list(12,"fourier"))
      b_1=fda::eval.basis(evalarg = seq(min(fts@grid[[j]][,1]),max(fts@grid[[j]][,1]),length.out=M_x),basisobj = fda::create.bspline.basis(rangeval = c(min(fts@grid[[j]][,1]),max(fts@grid[[j]][,1])),nbasis = B[[j]][[1]]))
      b_2=fda::eval.basis(evalarg = seq(min(fts@grid[[j]][,2]),max(fts@grid[[j]][,2]),length.out=M_y),basisobj = fda::create.fourier.basis(rangeval = c(min(fts@grid[[j]][,2]),max(fts@grid[[j]][,2])),nbasis = B[[j]][[2]]))
      fts@B[[j]]=kronecker(b_1,b_2)
      B[[j]]=fts@B[[j]]

    }else if(class(B[[j]])=="list" && B[[j]][[3]]=="fourier" && B[[j]][[4]]=="bspline" && length(grid[[j]])==2){
      # This case is similar to the previous else if statement
      # Example: list(list(10,"fourier"),list(12,"bspline"))
      b_1=fda::eval.basis(evalarg = seq(min(fts@grid[[j]][,1]),max(fts@grid[[j]][,1]),length.out=M_x),basisobj = fda::create.fourier.basis(rangeval = c(min(fts@grid[[j]][,1]),max(fts@grid[[j]][,1])),nbasis = B[[j]][[1]]))
      b_2=fda::eval.basis(evalarg = seq(min(fts@grid[[j]][,2]),max(fts@grid[[j]][,2]),length.out=M_y),basisobj = fda::create.bspline.basis(rangeval = c(min(fts@grid[[j]][,2]),max(fts@grid[[j]][,2])),nbasis = B[[j]][[2]]))
      fts@B[[j]]=kronecker(b_1,b_2)
      B[[j]]=fts@B[[j]]

    }else if(class(B[[j]])=="list" && B[[j]][[3]]=="fourier" && B[[j]][[4]]=="fourier" && length(grid[[j]])==2){
      # Generating a fourier basis for variables whose domain is two-dimensional using a supplied list.
      # Example: list(list(10,"fourier"),list(12,"fourier"))
      b_1=fda::eval.basis(evalarg = seq(min(fts@grid[[j]][,1]),max(fts@grid[[j]][,1]),length.out=M_x),basisobj = fda::create.fourier.basis(rangeval = c(min(fts@grid[[j]][,1]),max(fts@grid[[j]][,1])),nbasis = B[[j]][[1]]))
      b_2=fda::eval.basis(evalarg = seq(min(fts@grid[[j]][,2]),max(fts@grid[[j]][,2]),length.out=M_y),basisobj = fda::create.fourier.basis(rangeval = c(min(fts@grid[[j]][,2]),max(fts@grid[[j]][,2])),nbasis = B[[j]][[2]]))
      fts@B[[j]]=kronecker(b_1,b_2)
      B[[j]]=fts@B[[j]]

    }
    # Estimate the coefficients of each fts variable using basis elements and observed data.

    fts@C[[j]]=solve(t(B[[j]])%*%B[[j]])%*%t(B[[j]])%*%X[[j]]

    if(is.null(colnames(fts@C[[j]]))==TRUE){

      colnames(fts@C[[j]])=as.character(1:ncol(fts@C[[j]]))

    }

    if(is.null(rownames(fts@grid[[j]]))==TRUE){

      rownames(fts@grid[[j]])=as.character(1:nrow(fts@grid[[j]]))

    }

  }


  return(fts)

}
# class definitions, methods, etc. go outside of constructor.

check_fts <- function(object){

  errors <- character()
  msg <- NULL

  # check if fts lengths are equal

  if(length(object@C)>1){

    for(j_1 in 1:length(object@C)){
      for(j_2 in j_1:length(object@C)){


        if(length(dim(object@C[[j_1]]))==2 && length(dim(object@C[[j_2]]))==2){


          if(ncol(object@C[[j_1]])!=ncol(object@C[[j_2]])){

          msg<-paste("The length of variable", as.character(j_1),"and the length of variable", as.character(j_2),"are not equal.")

          errors <- c(errors,msg)

          msg <-NULL

          }


        }

        else if(length(dim(object@C[[j_1]]))==2 && length(dim(object@C[[j_2]]))==3){


          if(ncol(object@C[[j_1]])!=dim(object@C[[j_2]])[3]){

            msg<-paste("The length of variable", as.character(j_1),"and the length of variable", as.character(j_2),"are not equal.")

            errors <- c(errors,msg)

            msg <-NULL
          }


        }else if(length(dim(object@C[[j_1]]))==3 && length(dim(object@C[[j_2]]))==2){

          if(dim(object@C[[j_1]])[3]!=ncol(object@C[[j_2]])){

            msg<-paste("The length of variable", as.character(j_1),"and the length of variable", as.character(j_2),"are not equal.")

            errors <- c(errors,msg)

            msg <-NULL

          }

        }else if(length(dim(object@C[[j_1]]))==3 && length(dim(object@C[[j_2]]))==3){

          if(dim(object@C[[j_1]])[3]!=dim(object@C[[j_2]])[3]){

            msg<-paste("The length of variable", as.character(j_1),"and the length of variable", as.character(j_2),"are not equal.")

            errors <- c(errors,msg)

            msg <-NULL

          }

        }

      }
    }
  }

  # check if provided lists have same length
  if(length(object@C)!=length(object@B)){

    msg<-paste("The length of X is", as.character(length(object@C)),"and the length of B is", as.character(length(object@B)),"whereas the length of these lists should be equal.")

    errors <- c(errors,msg)

    msg <-NULL

  }
  if(length(object@C)!=length(object@grid)){

    msg<-paste("The length of X is", as.character(length(object@C)),"and the length of grid is", as.character(length(object@grid)),"whereas the length of these lists should be equal.")

    errors <- c(errors,msg)

    msg <-NULL

  }

  if(length(object@B)!=length(object@grid)){

    msg<-paste("The length of B is", as.character(length(object@B)),"and the length of grid is", as.character(length(object@grid)),"whereas the length of these lists should be equal.")

    errors <- c(errors,msg)

    msg <-NULL

  }
 # check if the classes and values for supplied inputs are correct
  for(i in 1L:length(object@C)){

    if (class(object@C[[i]])!="matrix" && class(object@C[[i]])!="numeric" && class(object@C[[i]])!="integer" && class(object@C[[i]])!="array"){

      msg <- paste("The class of entry", as.character(i) ,"of X is not matrix, numeric, or integer as is expected.")

      errors <- c(errors,msg)

      msg <-NULL
    }

  }

  for(i in 1L:length(object@B)){

    if (class(object@B[[i]])!="matrix" && class(object@B[[i]])!="numeric" && class(object@B[[i]])!="integer" && class(object@B[[i]])!="list"){

      msg <- paste("The class of entry", as.character(i) ,"of B is not matrix, numeric, integer, or list as is expected.")

      errors <- c(errors,msg)

      msg <-NULL

    }

    if(class(object@grid[[i]])=="list" && length(object@grid[[i]])==2){

      if (class(object@B[[i]])=="list" && length(object@B[[i]])!=4){

        msg <- paste("The length of entry", as.character(i), "of B is not four. The entry should be a basis matrix or a list of the form list(int_1,int_2,string_1,string_2) for variables taken over a two-dimensional domain where int_1,int_2>0 and string_1,string_2 = \"bspline\" or string_1,string_2 = \"fourier\".")

        errors <- c(errors,msg)

        msg <-NULL
      }

      if (class(object@B[[i]])=="list" && length(object@B[[i]])==4 && class(object@B[[i]][[1]])!="integer" && class(object@B[[i]][[1]])!="numeric"){


        msg <- paste("The first element in entry", as.character(i), "of B should be a positive integer.")

        errors <- c(errors,msg)

        msg <-NULL

      }

      if (class(object@B[[i]])=="list" && length(object@B[[i]])==4 && class(object@B[[i]][[1]])=="integer" || class(object@B[[i]][[1]])=="numeric" && object@B[[i]][[1]]<=0){

        msg <- paste("The first element in entry", as.character(i), "of B should be a positive integer.")

        errors <- c(errors,msg)

        msg <-NULL

      }


      if (class(object@B[[i]])=="list" && length(object@B[[i]])==4 && class(object@B[[i]][[2]])!="integer" && class(object@B[[i]][[2]])!="numeric"){

        msg <- paste("The second element in entry", as.character(i), "of B should be a positive integer.")

        errors <- c(errors,msg)

        msg <-NULL

      }

      if (class(object@B[[i]])=="list" && length(object@B[[i]])==4 && class(object@B[[i]][[2]])=="integer" || class(object@B[[i]][[2]])=="numeric" && object@B[[i]][[2]]<=0){

        msg <- paste("The second element in entry", as.character(i), "of B should be a positive integer.")

        errors <- c(errors,msg)

        msg <-NULL

      }


      if (class(object@B[[i]])=="list" && length(object@B[[i]])==4 && class(object@B[[i]][[3]])!="character"){

        msg <- paste("The third element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")

        errors <- c(errors,msg)

        msg <-NULL

      }

      if (class(object@B[[i]])=="list" && length(object@B[[i]])==4 && class(object@B[[i]][[3]])=="character" && object@B[[i]][[3]] != "bspline" && object@B[[i]][[3]]!= "fourier"){

        msg <- paste("The third element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")

        errors <- c(errors,msg)

        msg <-NULL

      }

      if (class(object@B[[i]])=="list" && length(object@B[[i]])==4 && class(object@B[[i]][[4]])!="character"){

        msg <- paste("The fourth element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")

        errors <- c(errors,msg)

        msg <-NULL

      }

      if (class(object@B[[i]])=="list" && length(object@B[[i]])==4 && class(object@B[[i]][[4]])=="character" && object@B[[i]][[4]] != "bspline" && object@B[[i]][[4]]!= "fourier"){

        msg <- paste("The fourth element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")

        errors <- c(errors,msg)

        msg <-NULL

      }




    }else{


      if (class(object@B[[i]])=="list" && length(object@B[[i]])!=2){

        msg <- paste("The length of entry", as.character(i), "of B is not two. The entry should be a basis matrix or a list of the form list(int,string) for variables taken over a one-dimensional domain where int>0 and string = \"bspline\" or string = \"fourier\".")

        errors <- c(errors,msg)

        msg <-NULL
      }

      if (class(object@B[[i]])=="list" && length(object@B[[i]])==2 && class(object@B[[i]][[1]])!="integer" && class(object@B[[i]][[1]])!="numeric"){

        msg <- paste("The first element in entry", as.character(i), "of B should be a positive integer.")

        errors <- c(errors,msg)

        msg <-NULL

      }

      if (class(object@B[[i]])=="list" && length(object@B[[i]])==2 && class(object@B[[i]][[1]])=="integer" || class(object@B[[i]][[1]])=="numeric" && object@B[[i]][[1]]<=0){

        msg <- paste("The first element in entry", as.character(i), "of B should be a positive integer.")

        errors <- c(errors,msg)

        msg <-NULL

      }

      if (class(object@B[[i]])=="list" && length(object@B[[i]])==2 && class(object@B[[i]][[2]])!="character"){

        msg <- paste("The second element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")

        errors <- c(errors,msg)

        msg <-NULL

      }

      if (class(object@B[[i]])=="list" && length(object@B[[i]])==2 && class(object@B[[i]][[2]])=="character" && object@B[[i]][[2]] != "bspline" && object@B[[i]][[2]]!= "fourier"){

        msg <- paste("The second element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")

        errors <- c(errors,msg)

        msg <-NULL

      }

    }

  }
 # Check classes and supplied inputs for the grids
  for(i in 1L:length(object@grid)){

    if(class(object@grid[[i]])=="list"){

      for(k in 1:length(object@grid[[i]])){

        if(class(object@grid[[i]][[k]])!="matrix" && class(object@grid[[i]][[k]])!="integer" && class(object@grid[[i]][[k]])!="numeric"){

          msg <- paste("Element ", as.character(k), " of entry ", as.character(i) ," of grid is not a matrix, numeric, or integer as is expected.",sep = "")
          errors <- c(errors,msg)

          msg <-NULL
        }

        if (class(object@grid[[i]][[k]])=="numeric" && object@grid[[i]][[k]][1]>=object@grid[[i]][[k]][2]){

          msg <- paste("The first element of interval ", as.character(k), " of entry ", as.character(i) ," of grid should be less than the second element of that same interval.",sep = "")

          errors <- c(errors,msg)

          msg <-NULL

        }


      }
    }else{

      if (class(object@grid[[i]])!="matrix" && class(object@grid[[i]])!="numeric" && class(object@grid[[i]])!="integer"){

        msg <- paste("The class of entry ", as.character(i) ," of grid is not matrix, numeric, integer, or list as is expected.",sep = "")

        errors <- c(errors,msg)

        msg <-NULL

      }


      if (class(object@grid[[i]])=="numeric" && object@grid[[i]][1]>=object@grid[[i]][2]){

        msg <- paste("The first element of entry ", as.character(i) ," of grid should be less than the second element.",sep = "")

        errors <- c(errors,msg)

        msg <-NULL

      }


    }



  }

  if(length(errors)!=0){
    return(errors)

  }
}
# Set the class definition
setClass("fts",representation(C = "list", B="list", grid="list"),validity=check_fts)
# Show method for the fts object
setMethod("show","fts",function(object){

  N=ncol(object@C[[1]])
  p=length(object@C)

  if(length(object@C)==1){

    cat("Univariate ", is(object), " of length ", N,".","\n",sep="")

  }else{

    cat("Multivariate ", is(object), " of length ", N, " that is comprised of ", p, " variables.","\n",sep="")


  }

  for(j in 1:p){

    if(ncol(object@grid[[j]])==1){
      domain_statement<-1
      object@grid[[j]]=matrix(data=object@grid[[j]],ncol=1)
      number_points<-nrow(object@grid[[j]])
      min_x=object@grid[[j]][1,1]
      max_x=object@grid[[j]][number_points,1]

      cat("\n","The domain of the observations for variable ", as.character(j), " is ", as.character(domain_statement), "-dimensional.","\n",

          "Observations for variable ", as.character(j) ," are sampled on a regular grid at ", as.character(number_points), " points.", "\n",

          "The x axis for this variable varies along the domain [",as.character(min_x),",",as.character(max_x),"].", "\n",




          sep = "")




    }else{

      domain_statement<-2
      number_points<-nrow(object@grid[[j]])
      min_x=object@grid[[j]][1,1]
      min_y=object@grid[[j]][1,2]
      max_x=object@grid[[j]][number_points,1]
      max_y=object@grid[[j]][number_points,2]

      cat("\n","The domain of the observations for variable ", as.character(j), " is ", as.character(domain_statement), "-dimensional.","\n",

          "Observations for variable ", as.character(j) ," are sampled on a regular grid at ", as.character(number_points), " points.", "\n",

          "The x axis for this variable varies along the domain [",as.character(min_x),",",as.character(max_x),"].", "\n",

          "The y axis for this variable varies along the domain [",as.character(min_y),",",as.character(max_y),"].", "\n",


          sep="")

    }



  }


})

