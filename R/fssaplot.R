#--------------------------------------------------------------
#' Plot Functional Singular Spectrum Analysis Objects
#'
#'  This is a plotting method for objects of class functional singular spectrum analysis (\code{\link{fssa}}). The method is designed to help the user make decisions
#'  on how to do the grouping stage of univariate or multivariate functional singular spectrum analysis.
#'
#' @param x an object of class \code{\link{fssa}}
#' @param d an integer which is the number of elementary components in the plot
#' @param idx a vector of indices of eigen elements to plot
#' @param idy a second vector of indices of eigen elements to plot (for type="paired")
#' @param groups a list or vector of indices determines grouping used for the decomposition(for type="wcor")
#' @param contrib a logical where if the value is 'TRUE' (the default), the contribution of the component to the total variance is displayed
#' @param type the type of plot to be displayed where possible types are:
#' \itemize{
#' \item \code{"values"} plot the square-root of singular values (default)
#' \item \code{"paired"} plot the pairs of eigenfunction's coefficients (useful for the detection of periodic components)
#' \item \code{"wcor"} plot the W-correlation matrix for the reconstructed objects
#' \item \code{"vectors"} plot the eigenfunction's coefficients (useful for the detection of period length)
#' \item \code{"lcurves"} plot of the eigenfunctions (useful for the detection of period length)
#' \item \code{"lheats"} heatmap plot the eigenfunctions (useful for the detection of meaningful patterns)
#' \item \code{"periodogram"} periodogram plot (useful for the detecting the frequencies of oscillations in functional data)
#' }
#' @param var an integer specifying the variable number
#' @param ylab the character vector of name of variables
#' @param ... arguments to be passed to methods, such as graphical parameters
#' @examples
#' \dontrun{
#' ## Simulated Data Example
#' require(Rfssa)
#' require(fda)
#' n <- 50 # Number of points in each function.
#' d <- 9
#' N <- 60
#' sigma <- 0.5
#' set.seed(110)
#' E <- matrix(rnorm(N*d,0,sigma/sqrt(d)),ncol = N, nrow = d)
#' basis <- create.fourier.basis(c(0, 1), d)
#' Eps <- fd(E,basis)
#' om1 <- 1/10
#' om2 <- 1/4
#' f0 <- function(tau, t) 2*exp(-tau*t/10)
#' f1 <- function(tau, t) 0.2*exp(-tau^3) * cos(2 * pi * t * om1)
#' f2 <- function(tau, t) -0.2*exp(-tau^2) * cos(2 * pi * t * om2)
#' tau <- seq(0, 1, length = n)
#' t <- 1:N
#' f0_mat <- outer(tau, t, FUN = f0)
#' f0_fd <- smooth.basis(tau, f0_mat, basis)$fd
#' f1_mat <- outer(tau, t, FUN = f1)
#' f1_fd <- smooth.basis(tau, f1_mat, basis)$fd
#' f2_mat <- outer(tau, t, FUN = f2)
#' f2_fd <- smooth.basis(tau, f2_mat, basis)$fd
#' Y_fd <- f0_fd+f1_fd+f2_fd
#' L <-10
#' U <- fssa(Y_fd,L)
#' plot(U)
#' plot(U,d=4,type="lcurves")
#' plot(U,d=4,type="vectors")
#' plot(U,d=5,type="paired")
#' plot(U,d=5,type="wcor")
#' plot(U,d=5,type="lheats")
#' plot(U,d=5,type="periodogram")
#' }
#' @seealso \code{\link{fssa}}, \code{\link{plot.fts}}
#' @note for a multivariate example, see the examples in \code{\link{fssa}}
#' @export
plot.fssa <- function(x, d = length(x$values),
                      idx = 1:d, idy = idx+1, contrib = TRUE,
                      groups = as.list(1:d), lagvec=1,
                      type = "values",var=1L,ylab=NA, animation = FALSE, ...) {



  rotate <- function(x) t(apply(x, 2, rev))

  p <- length(x$Y@C)
  A <- ((x$values)/sum(x$values))[1L:d]
  pr <- round(A * 100L, 2L)
  idx <- sort(idx)
  idy <- sort(idy)
  if(max(idx) > d | min(idx)< 1) stop("The idx must be subset of 1:d.")
  d_idx <- length(idx)
  if(contrib){
    main1 <- paste0(idx, "(", pr[idx],"%)")
    main2 <- paste0(idy, "(", pr[idy],"%)")
  } else{
    main1 <- paste(idx)
    main2 <- paste(idy)
  }
  N <- x$N
  L <- x$L
  K <- N-L+1L
  if (type %in% c("lheats","lcurves")) {
    if(ncol(x$Y@grid[[var]])==1){
      xindx <- x$Y@grid[[var]]
      z0 <- list()
      for (i in 1:d_idx){
        if(length(x$Y@C)==1){
          z0[[i]] <- x[[idx[i]]]
          }
        else{
          z0[[i]] <- x[[idx[i]]][[var]]
        }
      }
    }else{

      xindx <- unique(x$Y@grid[[var]][,"x"])
      yindx <- unique(x$Y@grid[[var]][,"y"])
      z0 <- list()
      for(i in 1:d_idx){
        z0_i<-list()
        for(j in 1:L){
          if(length(U$Y@C)==1){
            z0_i[[j]]<-matrix(data=x[[i]][,j],nrow=length(xindx),ncol=length(yindx))
          }else{
            z0_i[[j]]<-matrix(data=x[[i]][[var]][,j],nrow=length(xindx),ncol=length(yindx))
          }

        }
        z0[[i]]<-z0_i

      }

  }
  }
  if (type == "values") {


    val <- sqrt(x$values)[idx]
    graphics::plot(idx, val, type = "o", lwd = 3L,
         col = "dodgerblue3",
         ylab = "norms", xlab = "Components")
  } else if (type == "wcor") {
    W <- fwcor(x, groups)
    wplot(W)
  }  else  if (type == "lheats") {



    if(ncol(x$Y@grid[[var]])==1){
    n <- length(xindx)
    z <- c(sapply(z0, function(x) as.vector(x)))
    D0 <- expand.grid(x = 1L:L,
                      y = 1L:n, groups = idx)
    D0$z <- z
    D0$groups <- factor(rep(main1,
                               each = L * n), levels = main1)
    title0 <- "Singular functions"
    if(p>1) title0 <- paste(title0,"of the variable",
                            ifelse(is.na(ylab),var,ylab))
    p1 <- lattice::levelplot(z ~ x *
                               y | groups, data = D0,
                             xlab = "", ylab = "",
                             scales =  list(alternating=1,   # axes labels left/bottom
                                            tck = c(1,0)), as.table = TRUE,
                             main = title0,
                             col.regions = grDevices::heat.colors(100))
    graphics::plot(p1)
    }else{

      if(animation==FALSE){

      par(ask=TRUE)
      for(j in 1:L){

        print(lattice::levelplot(z0[[lagvec]][[j]]~x$Y@grid[[var]][,"x"]*x$Y@grid[[var]][,"y"],xlab="",ylab="",col.regions = heat.colors(100),
                           main=paste("Function ",as.character(j)," of lag vector ",as.character(lagvec),sep="")))
      }

      }else{


        if(length(x$Y@C)!=1){
          # MFSSA animation images
          tickfont <- 1.7 # font size of colorbar tick marks
          titlefont <- 2.5 # font size of main title
          subtitlefont <- 2 # font size of sub plots
          axisfont<-1.9
          labelfont<-2
          par(ask=TRUE)
          pdf(file="./mfssa_plots/mfssa_NDVI_fun_ani.pdf",width=6,height=6);
        for(j in 1:L){
          d1 <- floor(sqrt(d_idx))
          d2 <- ceiling(d_idx/d1)
          graphics::par(mfrow = c(d1, d2),
                        mar = c(0, 1, 7, 2),oma=c(5,5,3,1),cex.main=titlefont)
          title0 <- paste("(B) Function " ,as.character(j)," of Lag Vectors",sep="")
          for (i in 1:d_idx){
            if(i==1){
            graphics::image(x = 1:33,y=1:33,z=rotate(matrix(data=x[[i]][[var]][,j],nrow=33,ncol=33)),useRaster=TRUE,
                            xlab="lon",ylab="lat",main=main1[i], xaxt="n", yaxt="n",
                            cex.axis = axisfont, cex.main=titlefont, cex.lab=labelfont,
                            )
              #axis(side = 1, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.303323,113.407503,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))
              axis(side = 2, las=2,cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(48.706792,48.775542,length.out=4),decreasing = FALSE)),start=1,stop=5),at = seq(1,33,length.out=4))
              axis(side = 3, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.458407,113.56273,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))



            }else if(i==2){


            graphics::image(x = 1:33,y=1:33,z=rotate(matrix(data=x[[i]][[var]][,j],nrow=33,ncol=33)),useRaster=TRUE,
                            xlab="Latitude",ylab="Longitude",main=main1[i], xaxt="n", yaxt="n",
                            cex.axis = axisfont, cex.main=titlefont, cex.lab=labelfont,
            )

              #axis(side = 1, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.303323,113.407503,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))
              #axis(side = 2, las=2,cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(48.706792,48.775542,length.out=4),decreasing = FALSE)),start=1,stop=5),at = seq(1,33,length.out=4))
              axis(side = 3, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.458407,113.56273,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))





          }else if(i==3){
            graphics::image(x = 1:33,y=1:33,z=rotate(matrix(data=x[[i]][[var]][,j],nrow=33,ncol=33)),useRaster=TRUE,
                            xlab="lat",ylab="lon",main=main1[i], xaxt="n", yaxt="n",
                            cex.axis = axisfont, cex.main=titlefont, cex.lab=labelfont,
            )

            axis(side = 1, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.303323,113.407503,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))
            axis(side = 2, las=2,cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(48.706792,48.775542,length.out=4),decreasing = FALSE)),start=1,stop=5),at = seq(1,33,length.out=4))
            #axis(side = 3, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.458407,113.56273,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))




          }else if(i==4){
            graphics::image(x = 1:33,y=1:33,z=rotate(matrix(data=x[[i]][[var]][,j],nrow=33,ncol=33)),useRaster=TRUE,
                            xlab="Latitude",ylab="Longitude",main=main1[i],  xaxt="n", yaxt="n",
                            cex.axis = axisfont, cex.main=titlefont, cex.lab=labelfont,
            )
            graphics::title(title0,outer = TRUE,xlab="lon",ylab="lat",cex=titlefont,cex.lab=labelfont)

            axis(side = 1, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.303323,113.407503,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))
            #axis(side = 2, las=2,cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(48.706792,48.775542,length.out=4),decreasing = FALSE)),start=1,stop=5),at = seq(1,33,length.out=4))
            #axis(side = 3, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.458407,113.56273,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))



          }
          }
        }
          dev.off()
        }else{
          # UFSSA on images animation (A)
          tickfont <- 1.7 # font size of colorbar tick marks
          titlefont <- 2.5 # font size of main title
          subtitlefont <- 2 # font size of sub plots
          axisfont<-1.9
          labelfont<-2
          par(ask=TRUE)
          pdf(file="./mfssa_plots/NDVI_fun_ani.pdf",width=6,height=6);
          for(j in 1:L){
            d1 <- floor(sqrt(d_idx))
            d2 <- ceiling(d_idx/d1)
            graphics::par(mfrow = c(d1, d2),
                          mar = c(0, 1, 7, 2),oma=c(5,5,3,1),cex.main=titlefont)
            title0 <- paste("(A) Function " ,as.character(j)," of Lag Vectors",sep="")
            for (i in 1:d_idx){
              if(i==1){
                graphics::image(x = 1:33,y=1:33,z=rotate(matrix(data=x[[i]][,j],nrow=33,ncol=33)),useRaster=TRUE,
                                xlab="Latitude",ylab="Longitude",main=main1[i], xaxt="n",yaxt="n",
                                cex.axis = axisfont, cex.main=titlefont, cex.lab=labelfont,
                )
                #axis(side = 1, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.303323,113.407503,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))
                axis(side = 2, las=2,cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(48.706792,48.775542,length.out=4),decreasing = FALSE)),start=1,stop=5),at = seq(1,33,length.out=4))
                axis(side = 3, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.458407,113.56273,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))


              }else if(i==2){


                graphics::image(x = 1:33,y=1:33,z=rotate(matrix(data=x[[i]][,j],nrow=33,ncol=33)),useRaster=TRUE,
                                xlab="Latitude",ylab="Longitude",main=main1[i], xaxt="n", yaxt="n",
                                cex.axis = axisfont, cex.main=titlefont, cex.lab=labelfont,
                )

                #axis(side = 1, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.303323,113.407503,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))
                #axis(side = 2, las=2,cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(48.706792,48.775542,length.out=4),decreasing = FALSE)),start=1,stop=5),at = seq(1,33,length.out=4))
                axis(side = 3, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.458407,113.56273,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))




              }else if(i==3){
                graphics::image(x = 1:33,y=1:33,z=rotate(matrix(data=x[[i]][,j],nrow=33,ncol=33)),useRaster=TRUE,
                                xlab="Latitude",ylab="Longitude",main=main1[i], xaxt="n",yaxt="n",
                                cex.axis = axisfont, cex.main=titlefont, cex.lab=labelfont,
                )

                axis(side = 1, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.303323,113.407503,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))
                axis(side = 2, las=2,cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(48.706792,48.775542,length.out=4),decreasing = FALSE)),start=1,stop=5),at = seq(1,33,length.out=4))
                #axis(side = 3, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.458407,113.56273,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))




              }else if(i==4){
                graphics::image(x = 1:33,y=1:33,z=rotate(matrix(data=x[[i]][,j],nrow=33,ncol=33)),useRaster=TRUE,
                                xlab="Latitude",ylab="Longitude",main=main1[i], xaxt="n",yaxt="n",
                                cex.axis = axisfont, cex.main=titlefont, cex.lab=labelfont,
                )

                axis(side = 1, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.303323,113.407503,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))
                #axis(side = 2, las=2,cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(48.706792,48.775542,length.out=4),decreasing = FALSE)),start=1,stop=5),at = seq(1,33,length.out=4))
                #axis(side = 3, cex.lab=flab,cex.axis=ftick,labels = substr(x=as.character(sort(x=seq(113.458407,113.56273,length.out=4),decreasing = TRUE)),start=1,stop=6),at = seq(1,33,length.out=4))
                graphics::title(title0,outer = TRUE,xlab="lon",ylab="lat",cex=titlefont,cex.lab=labelfont)
              }
            }

          }
          dev.off()




        }



      }

    }
  } else if (type == "lcurves") {
    tickfont <- 1.6 # font size of colorbar tick marks
    titlefont <- 2.2 # font size of main title
    subtitlefont <- 1.9 # font size of sub plots
    axisfont<-1.9
    labelfont<-2
    if(ncol(x$Y@grid[[var]])==1){
      if(animation==FALSE){
    col2 <- grDevices::rainbow(L)
    d1 <- floor(sqrt(d_idx))
    d2 <- ceiling(d_idx/d1)
    pdf(file="./mfssa_plots/temp_funs.pdf",width=6,height=6);
    graphics::par(mfrow = c(d1, d2),
                  mar = c(0.25, .5, 2, 0),oma=c(5,5,3,1),cex.main=titlefont)
    title0 <- "(A) Singular Functions"
    if(p>1) title0 <- paste(title0,"of the variable",
                            ifelse(is.na(ylab),var,ylab))
    if(length(x$Y@C)==1){

      for (i in 1:d_idx){
        if(i==1){
          graphics::matplot(x[[idx[i]]],
                         lty = 1, xlab = "",ylim=range(z0),
                         main = main1[i], ylab = "",
                         lwd = 3, col = col2,xaxt="n",cex.axis=tickfont,type="l")
        }else if(i==2){
          #par(mar=c(5.1,4.1,3,2.1))
          graphics::matplot(x[[idx[i]]],
                         lty = 1, xlab = "",ylim=range(z0),
                         main = main1[i], ylab = "",
                         lwd = 3, col = col2,xaxt="n",yaxt="n",cex.axis=tickfont,type="l")


        }else if(i==3){
          graphics::matplot(x[[idx[i]]],
                         lty = 1, xlab = "",ylim=range(z0),
                         main = main1[i], ylab = "",
                         lwd = 3, col = col2,cex.axis=tickfont,type="l")


        }else{
          graphics::matplot(x[[idx[i]]],
                         lty = 1, xlab = "",ylim=range(z0),
                         main = main1[i], ylab = "",
                         lwd = 3, col = col2,yaxt="n",cex.axis=tickfont,type="l")


          graphics::title(main=title0,xlab="Intraday Time",ylab="Temperature (C)",cex.lab=labelfont,outer=TRUE)
        }

      }
      dev.off()

    }else{
      pdf(file="./mfssa_plots/mfssa_temp_funs.pdf",width=6,height=6);
      graphics::par(mfrow = c(d1, d2),
                    mar = c(1, .5, 2.5, 0),oma=c(5,5,3,1),cex.main=titlefont)
      title0 <- "(B) Singular Functions"
      for (i in 1:d_idx){
        if(i==1){
          graphics::matplot(x[[idx[i]]][[var]],
                         lty = 1, xlab = "",ylim=range(z0),
                         main = main1[i], ylab = "",
                         lwd = 3, col = col2,xaxt="n",cex.axis=tickfont,type="l")
        }else if(i==2){
          #par(mar=c(5.1,4.1,3,2.1))
          graphics::matplot(x[[idx[i]]][[var]],
                         lty = 1, xlab = "",ylim=range(z0),
                         main = main1[i], ylab = "",
                         lwd = 3, col = col2,xaxt="n",yaxt="n",cex.axis=tickfont,type="l")



        }else if(i==3){
          graphics::matplot(x[[idx[i]]][[var]],
                         lty = 1, xlab = "",ylim=range(z0),
                         main = main1[i], ylab = "",
                         lwd = 3, col = col2,cex.axis=tickfont,type="l")


        }else{
          graphics::matplot(x[[idx[i]]][[var]],
                         lty = 1, xlab = "",ylim=range(z0),
                         main = main1[i], ylab = "",
                         lwd = 3, col = col2,yaxt="n",cex.axis=tickfont,type="l")




          graphics::title(main=title0,xlab="Intraday Time",ylab="Normalized Temperature",outer=TRUE,cex.lab=labelfont)
        }

      }
      dev.off()

    }

      }else{
        # MFSSA curves
        tickfont <- 1.7 # font size of colorbar tick marks
        titlefont <- 2.2 # font size of main title
        subtitlefont <- 1.9 # font size of sub plots
        axisfont<-1.9
        labelfont<-2
        if(length(x$Y@C)!=1){
        par(ask=TRUE)
        pdf(file="./mfssa_plots/mfssa_temp_funs_ani.pdf",width=6,height=6);
        for (j in 1:L){
          col2 <- grDevices::rainbow(L)
          d1 <- floor(sqrt(d_idx))
          d2 <- ceiling(d_idx/d1)
          graphics::par(mfrow = c(d1, d2),
                        mar = c(1, .5, 2.5, 0),oma=c(5,5,3,1),cex.main=titlefont)
          title0 <- paste("(D) Function ",as.character(j)," of Lag Vectors",sep="")
          for(i in 1:d_idx){
            if(i==1){
              graphics::matplot(x[[idx[i]]][[var]][,j],
                                lty = 1, xlab = "",ylim=range(z0),
                                main = main1[i], ylab = "",
                                lwd = 3, col = col2,xaxt="n",cex.axis=tickfont,type="l")
            }else if(i==2){
              #par(mar=c(5.1,4.1,3,2.1))
              graphics::matplot(x[[idx[i]]][[var]][,j],
                                lty = 1, xlab = "",ylim=range(z0),
                                main = main1[i], ylab = "",
                                lwd = 3, col = col2,xaxt="n",yaxt="n",cex.axis=tickfont,type="l")


            }else if(i==3){
              graphics::matplot(x[[idx[i]]][[var]][,j],
                                lty = 1, xlab = "",ylim=range(z0),
                                main = main1[i], ylab = "",
                                lwd = 3, col = col2,cex.axis=tickfont,type="l")


            }else{
              graphics::matplot(x[[idx[i]]][[var]][,j],
                                lty = 1, xlab = "",ylim=range(z0),
                                main = main1[i], ylab = "",
                                lwd = 3, col = col2,yaxt="n",cex.axis=tickfont,type="l")


              graphics::title(main=title0,xlab="Intraday Time",ylab="Normalized Temperature",cex.lab=labelfont,outer=TRUE)
            }
          }
        }
        dev.off()
        }else{

          #par(ask=TRUE) FSSA temp animation
          # here's what I want for animation
          pdf(file="./mfssa_plots/temp_funs_ani.pdf",width=6,height=6);
          for (j in 1:L){
            col2 <- grDevices::rainbow(L)
            d1 <- floor(sqrt(d_idx))
            d2 <- ceiling(d_idx/d1)
            graphics::par(mfrow = c(d1, d2),
                          mar = c(0.25, .5, 2, 0),oma=c(5,5,3,1),cex.main=titlefont)
            title0 <- paste("(C) Function ",as.character(j)," of Lag Vectors",sep="")
            for(i in 1:d_idx){
              if(i==1){
                graphics::matplot(x[[idx[i]]][,j],
                                  lty = 1, xlab = "",ylim=range(z0),
                                  main = main1[i], ylab = "",
                                  lwd = 3, col = col2,xaxt="n",cex.axis=tickfont,type="l")
              }else if(i==2){
                #par(mar=c(5.1,4.1,3,2.1))
                graphics::matplot(x[[idx[i]]][,j],
                                  lty = 1, xlab = "",ylim=range(z0),
                                  main = main1[i], ylab = "",
                                  lwd = 3, col = col2,xaxt="n",yaxt="n",cex.axis=tickfont,type="l")


              }else if(i==3){
                graphics::matplot(x[[idx[i]]][,j],
                                  lty = 1, xlab = "",ylim=range(z0),
                                  main = main1[i], ylab = "",
                                  lwd = 3, col = col2,cex.axis=tickfont,type="l")


              }else{
                graphics::matplot(x[[idx[i]]][,j],
                                  lty = 1, xlab = "",ylim=range(z0),
                                  main = main1[i], ylab = "",
                                  lwd = 3, col = col2,yaxt="n",cex.axis=tickfont,type="l")


                graphics::title(main=title0,xlab="Intraday Time",ylab="Temperature (C)",cex.lab=labelfont,outer=TRUE)
              }
            }
          }
          dev.off()


        }


      }

    }else{

      stop("lcurves defined only for variables whose observations are curves")

    }
    graphics::par(mfrow = c(1, 1))
  } else if (type == "vectors"){
    titlefont <- 2 # font size of main title
    subtitlefont <- 1.9 # font size of sub plots
    x0 <- c(apply(x$RVectrs[,idx],2,scale,center=F))
    D0 <- data.frame(x = x0,
                     time = rep(1L:K, d_idx))
    D0$groups <- factor(rep(main1,
                               each = K), levels = main1)
    p1 <- lattice::xyplot(x ~ time |
                            groups, data = D0, par.strip.text=list(cex=subtitlefont), xlab = "",
                          ylab = "", main = list("(C) Singular Vectors",cex=2),lwd=3,
                          scales = list(x = list(cex=c(1.4,1.4)),
                                        y = list(cex=c(1.4, 1.4), # increase font size
                                                 alternating=1,   # axes labels left/bottom
                                                 tck = c(1,0))),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
  } else if (type == "paired"){
    titlefont <- 1.3 # font size of main title
    subtitlefont <- 1.2 # font size of sub plots
    axisfont<-1.9
    d_idy <- length(idy)
    if(d_idx != d_idy) stop("The length of idx and idy must be same")
    x0 <- c(apply(x$RVectrs[,idx],2,scale,center=F))
    y0 <- c(apply(x$RVectrs[,idy],2,scale,center=F))
    D0 <- data.frame(x = x0, y = y0)
    main3 <- paste(main1, "vs", main2)
    D0$groups <- factor(rep(main3, each = K), levels = main3)
    p1 <- lattice::xyplot(x ~ y | groups,
                          data = D0, xlab = "", par.strip.text=list(cex=subtitlefont),
                          ylab = "", main = "Paired Singular vectors (Right)",
                          scales = list(x = list(at = NULL, relation="same"),
                                        y = list(at = NULL, relation="same")),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
  } else if (type == "periodogram"){
    ff <- function(x) {
      I <- abs(fft(x)/sqrt(K))^2
      P = (4/K) * I
      return(P[1:(floor(K/2) + 1)])
    }
    x0 <- c(apply(apply(x$RVectrs[,idx],2,scale,center=F),2,ff))
    D0 <- data.frame(x = x0,
                     time = rep((0:floor(K/2))/K, d_idx))
    D0$groups <- factor(rep(main1,
                               each = (floor(K/2) + 1)), levels = main1)
    p1 <- lattice::xyplot(x ~ time |
                            groups, data = D0, xlab = "",
                          ylab = "", main = "Periodogram of Singular vectors",
                          scales = list(y = list(at = NULL, relation="same")),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
    }else {
    stop("Unsupported type of fssa plot!")
    }



}
