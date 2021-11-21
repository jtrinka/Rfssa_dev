#--------------------------------------------------------------
#' Functional Time Series Visualization Tools Using Plotly
#'
#' This is a plotting method for univariate or multivariate functional time series (\code{\link{fts}}). This method is designed to help the user visualize
#' \code{\link{fts}} data using a variety of techniques that use plotly.
#'
#' @param x an object of class \code{\link{fts}}
#' @param types the types of plot to be displayed where possible types are:
#' \itemize{
#' \item \code{"line"} plot the \code{\link{fts}} elements in a line plot (default)
#' \item \code{"heatmap"} plot the \code{\link{fts}} elements in a heat map
#' \item \code{"3Dsurface"} plot the \code{\link{fts}} elements as a surface
#' \item \code{"3Dline"} plot the \code{\link{fts}} elements in a three-dimensional line plot
#' }
#' @param npts number of points to evaluate functional object at
#' @param main the main title
#' @param ylab the y-axis label
#' @param xlab the x-axis label
#' @param tlab the time-axis label
#' @param var an integer specifying the variable number to plot if \code{types="3Dsurface"} or \code{types="3Dline"}
#' @param ... arguments to be passed to methods, such as graphical parameters.
#' @importFrom plotly plot_ly add_lines layout subplot add_surface hide_colorbar
#' @importFrom hrbrthemes theme_ipsum
#' @import dplyr
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
#' plot(Y,types = "heatmap")
#' plot(Y,types = "line",var = 1)
#' plot(Y,types = "3Dsurface",var = 1)
#' plot(Y,types = "3Dline",var = 1)
#' }
#'
#' @note for a multivariate example, see the examples in \code{\link{fssa}}
#'
#'
#' @export
plot.fts <- function(Y,vars=NULL,types=NULL,subplot=TRUE,main=NULL,ylabels=NULL,xlabels=NULL,tlabels=NULL,zlabels=NULL, ...){

  p <- length(Y@C)
  N <- ncol(Y@C[[1]])
  time <- colnames(Y@C[[1]])
  count_twod <- 0
  Pl <- list()

  if(is.null(types)) types=rep(NA,p);
  if(is.null(ylabels)) ylabels=rep(NA,p);
  if(is.null(xlabels)) xlabels=rep(NA,p);
  if(is.null(tlabels)) tlabels=rep(NA,p);
  if(is.null(zlabels)) zlabels=rep(NA,p);
  if(is.null(main)) main=rep(NA,p);
  if(is.null(vars) == FALSE && length(types)!=length(vars)) return(stop("\"vars\" and \"types\" should be the same length."));


  if(is.null(vars)==TRUE){
    for(j in 1:p){

      if(ncol(Y@grid[[j]])==1){

        if(is.na(types[j])==TRUE||types[j]=="line"){

          if(is.na(ylabels[j])) ylabels[j] <- "y";
          if(is.na(xlabels[j])) xlabels[j] <- "x";
          if(is.na(main[j]) && p==1 || is.na(main[j]) && subplot==FALSE) main[j] <- paste("Variable",j);
          if(subplot==TRUE && length(types)>1) main[j]=NA
          y <- tibble::as_tibble(data.frame(y=c(Y@B[[j]]%*%Y@C[[j]])))
          y$time <- as.factor(rep(time,each=nrow(Y@grid[[j]])))
          y$x <- rep(1:nrow(Y@grid[[j]]),ncol(Y@C[[j]]))
          Pl[[j]] <- y %>%
          group_by(time) %>%
          plot_ly(x=~x,y=~y) %>%
          add_lines(color = ~time,colors=c("lightsteelblue","royalblue4"),showlegend=FALSE) %>%
          layout(title=main[j],xaxis = list(title = xlabels[j]),yaxis = list(title =ylabels[j]))


        }else if(types[j]=="heatmap"){


          if(is.na(ylabels[j])) ylabels[j] <- "x";
          if(is.na(xlabels[j])) tlabels[j] <- "t";
          if(is.na(main[j]) && p==1 || is.na(main[j]) && subplot==FALSE) main[j] <- paste("Variable",j);
          if(subplot==TRUE && length(types)>1) main[j]=NA
          z0 <- Y@B[[j]]%*%Y@C[[j]]
          Pl[[j]] <- plot_ly(z = z0, x=time, y = u, types = "heatmap", colorscale = list(c(0,'#FFFFFAFF'), c(1,'#FF0000FF')),
                             showscale =FALSE) %>%
            layout(title = main[j], yaxis = list(title = ylabels[j]),xaxis = list(title = tlabels[j]))
        }else if(types[j]=="3Dsurface"){

          z0 <- Y@B[[j]]%*%Y@C[[j]]
          axx <-axy <-axz <- list(
            gridcolor="rgb(180, 180, 180)",
            zerolinecolor="rgb(255,255,255)"
          )
          axx$title <- ifelse(is.na(tlabels[j]),"t",tlabels[j])
          axy$title <- ifelse(is.na(xlabels[j]),"x",xlabels[j])
          axz$title <- ifelse(is.na(ylabels[j]),paste("Variable",j),ylabels[j])
         Pl[[j]] <- plot_ly(z = z0, x=time, y = u, colorscale = list(c(0,'#FFFFFAFF'), c(1,'#FF0000FF'))) %>%
            layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))%>%
            add_surface(showscale=FALSE)

        }else if (types[j]=="3Dline"){

          D0 <- tibble::as_tibble(data.frame(z=c(Y@B[[j]]%*%Y@C[[j]])))
          D0$time <- rep(time,each=nrow(Y@grid[[j]]))
          D0$x <- rep(1:nrow(Y@grid[[j]]),ncol(Y@C[[j]]))
          axx <-axy <-axz <- list(
            gridcolor="rgb(180, 180, 180)",
            zerolinecolor="rgb(255,255,255)"
          )
          axx$title <- ifelse(is.na(tlabels[j]),"t",tlabels[j])
          axy$title <- ifelse(is.na(xlabels[j]),"x",xlabels[j])
          axz$title <- ifelse(is.na(ylabels[j]),paste("Variable",j),ylabels[j])
          Pl[[j]]<- D0 %>%
            group_by(time) %>%
            plot_ly(x=~time,z=~z,y=~x, type = 'scatter3d', mode = 'lines', color = ~z,
                    line = list(width = 4), colors=c("#FFFFFAFF","#FF0000FF")) %>%
            layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz)) %>% hide_colorbar()


        }

      }else{
        if(is.na(types[j])==FALSE && types[j]!="heatmap") cat("Notice: The only plotting option available for variables observed over two-dimensional domains is \"heatmap\". Future plotting options will be added for these types of variables in the future.");
        count_twod=count_twod+1
        if(is.na(ylabels[j])) ylabels[j] <- "y";
        if(is.na(xlabels[j])) xlabels[j] <- "x";
        if(is.na(zlabels[j])) zlabels[j] <- "z";
        if(is.na(main[j])) main[j] <- paste("Variable",j);
        y <- tibble::as_tibble(data.frame(y=c(Y@B[[j]]%*%Y@C[[j]])))
        time=as.character(1:ncol(Y@C[[1]]))
        y$time <- as.factor(rep(time,each=nrow(Y@grid[[j]])))
        y$x_1 <- rep(Y@grid[[j]][,1],ncol(Y@C[[j]]))
        y$x_2 <- rep(Y@grid[[j]][,2],ncol(Y@C[[j]]))
        Pl[[j]] <- ggplotly(ggplot(y,aes(x_1,x_2,fill=y,frame=time))+geom_tile()+scale_fill_distiller(palette = "RdYlBu")+theme_ipsum()+xlab(xlabels[j])+
                              ylab(ylabels[j])+labs(fill=zlabels[j])+ggtitle(main[j]))
        }


      }

  }else{


    for(j in 1:length(vars)){

      if(ncol(Y@grid[[vars[j]]])==1){

        if(is.na(types[j])==TRUE||types[j]=="line"){

          if(is.na(ylabels[j])) ylabels[j] <- "y";
          if(is.na(xlabels[j])) xlabels[j] <- "x";
          if(is.na(main[j]) && length(vars)==1 || is.na(main[j]) && subplot==FALSE) main[j] <- paste("Variable",vars[j]);
          if(subplot==TRUE && length(vars)>1) main[j]=NA
          y <- tibble::as_tibble(data.frame(y=c(Y@B[[vars[j]]]%*%Y@C[[vars[j]]])))
          y$time <- as.factor(rep(time,each=nrow(Y@grid[[vars[j]]])))
          y$x <- rep(1:nrow(Y@grid[[vars[j]]]),ncol(Y@C[[vars[j]]]))
          Pl[[j]] <- y %>%
            group_by(time) %>%
            plot_ly(x=~x,y=~y) %>%
            add_lines(color = ~time,colors=c("lightsteelblue","royalblue4"),showlegend=FALSE) %>%
            layout(title=main[j],yaxis = list(title = ylabels[j]),xaxis = list(title = xlabels[j]))


        }else if(types[j]=="heatmap"){

          if(is.na(ylabels[j])) ylabels[j] <- "x";
          if(is.na(xlabels[j])) tlabels[j] <- "t";
          if(is.na(main[j]) && length(vars)==1 || is.na(main[j]) && subplot==FALSE) main[j] <- paste("Variable",vars[j]);
          if(subplot==TRUE && length(vars)>1) main[j]=NA
          z0 <- Y@B[[vars[j]]]%*%Y@C[[vars[j]]]
          Pl[[j]] <- plot_ly(z = z0, x=time, y = u, types = "heatmap", colorscale = list(c(0,'#FFFFFAFF'), c(1,'#FF0000FF')),
                             showscale =FALSE) %>%
            layout(title = main[j], yaxis = list(title = ylabels[j]),xaxis = list(title = tlabels[j]))
        }else if(types[j]=="3Dsurface"){

          z0 <- Y@B[[vars[j]]]%*%Y@C[[vars[j]]]
          axx <-axy <-axz <- list(
            gridcolor="rgb(180, 180, 180)",
            zerolinecolor="rgb(255,255,255)"
          )
          axx$title <- ifelse(is.na(tlabels[j]),"t",tlabels[j])
          axy$title <- ifelse(is.na(xlabels[j]),"x",xlabels[j])
          axz$title <- ifelse(is.na(ylabels[j]),paste("Variable",vars[j]),ylabels[j])
          Pl[[j]] <- plot_ly(z = z0, x=time, y = u, colorscale = list(c(0,'#FFFFFAFF'), c(1,'#FF0000FF'))) %>%
            layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))%>%
            add_surface(showscale=FALSE)

        }else if (types[j]=="3Dline"){

          D0 <- tibble::as_tibble(data.frame(z=c(Y@B[[vars[j]]]%*%Y@C[[vars[j]]])))
          D0$time <- rep(time,each=nrow(Y@grid[[vars[j]]]))
          D0$x <- rep(1:nrow(Y@grid[[vars[j]]]),ncol(Y@C[[vars[j]]]))
          axx <-axy <-axz <- list(
            gridcolor="rgb(180, 180, 180)",
            zerolinecolor="rgb(255,255,255)"
          )
          axx$title <- ifelse(is.na(tlabels[j]),"t",tlabels[j])
          axy$title <- ifelse(is.na(xlabels[j]),"x",xlabels[j])
          axz$title <- ifelse(is.na(ylabels[j]),paste("Variable",vars[j]),ylabels[j])
          Pl[[j]]<- D0 %>%
            group_by(time) %>%
            plot_ly(x=~time,z=~z,y=~x, type = 'scatter3d', mode = 'lines', color = ~z,
                    line = list(width = 4), colors=c("#FFFFFAFF","#FF0000FF")) %>%
            layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz)) %>% hide_colorbar()


        }

      }else{
        if(is.na(types[j])==FALSE && types[j]!="heatmap") cat("Notice: The only plotting option available for variables observed over two-dimensional domains is \"heatmap\". Future plotting options will be added for these types of variables in the future.");
        count_twod=count_twod+1
        if(is.na(ylabels[j])) ylabels[j] <- "y";
        if(is.na(xlabels[j])) xlabels[j] <- "x";
        if(is.na(zlabels[j])) zlabels[j] <- "z";
        if(is.na(main[j])) main[j] <- paste("Variable",vars[j]);
        y <- tibble::as_tibble(data.frame(y=c(Y@B[[vars[j]]]%*%Y@C[[vars[j]]])))
        time=as.character(1:ncol(Y@C[[1]]))
        y$time <- as.factor(rep(time,each=nrow(Y@grid[[vars[j]]])))
        y$x_1 <- rep(Y@grid[[vars[j]]][,1],ncol(Y@C[[vars[j]]]))
        y$x_2 <- rep(Y@grid[[vars[j]]][,2],ncol(Y@C[[vars[j]]]))
        Pl[[j]] <- ggplotly(ggplot(y,aes(x_1,x_2,fill=y,frame=time))+geom_tile()+scale_fill_distiller(palette = "RdYlBu")+theme_ipsum()+xlab(xlabels[j])+
                              ylab(ylabels[j])+labs(fill=zlabels[j])+ggtitle(main[j]))

      }


    }



  }

  if(count_twod>=1){

    print(Pl)

  }else{

    if("3Dsurface"%in%types || "3Dline"%in%types){

      print(Pl)

    }else if("3Dsurface"%in%types==FALSE && subplot!=TRUE || "3Dline"%in%types==FALSE && subplot!=TRUE){

      print(Pl)

    }else{


      print(subplot(Pl))

    }
  }
}
