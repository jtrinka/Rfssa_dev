#--------------------------------------------------------------
#' Functional Time Series Visualization Tools Using Plotly
#'
#' This is a plotting method for univariate or multivariate functional time series (\code{\link{fts}}). This method is designed to help the user visualize
#' \code{\link{fts}} data using a variety of techniques that use plotly.
#'
#' @param x an object of class \code{\link{fts}}
#' @param type the type of plot to be displayed where possible types are:
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
#' @param var an integer specifying the variable number to plot if \code{type="3Dsurface"} or \code{type="3Dline"}
#' @param ... arguments to be passed to methods, such as graphical parameters.
#' @importFrom plotly plot_ly add_lines layout subplot add_surface hide_colorbar
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
#' plot(Y,type = "heatmap")
#' plot(Y,type = "line",var = 1)
#' plot(Y,type = "3Dsurface",var = 1)
#' plot(Y,type = "3Dline",var = 1)
#' }
#'
#' @note for a multivariate example, see the examples in \code{\link{fssa}}
#'
#'
#' @export
plot.fts <- function(Y,vars=NULL,type=NULL,subplot=TRUE,main=NULL,ylab="y",xlab="x",tlab="t",zlab="z", ...){

  p <- length(Y@C)
  N <- ncol(Y@C[[1]])
  time <- colnames(Y@C[[1]])
  Pl <- list()

  if(is.null(type)==TRUE){

    type=rep(NA,p)

  }

  if(is.null(vars)==TRUE){
    for(j in 1:p){

      if(ncol(Y@grid[[j]])==1){

        if(is.na(type[j])==TRUE||type[j]=="line"){

          y <- tibble::as_tibble(data.frame(y=c(Y@B[[j]]%*%Y@C[[j]])))
          y$time <- as.factor(rep(time,each=nrow(Y@grid[[j]])))
          y$x <- rep(1:nrow(Y@grid[[j]]),ncol(Y@C[[j]]))
          if(ylab=="y") y_var <- paste("Variable",j) else y_var <- ylab[j];
          Pl[[j]] <- y %>%
          group_by(time) %>%
          plot_ly(x=~x,y=~y) %>%
          add_lines(color = ~time,colors=c("lightsteelblue","royalblue4"),showlegend=FALSE) %>%
          layout(title=main,yaxis = list(title = ylab),xaxis = list(title = xlab))


        }else if(type[j]=="heatmap"){

          z0 <- Y@B[[j]]%*%Y@C[[j]]
          if(ylab=="y") y_var <- paste("Variable",j) else y_var <- ylab[j];
          Pl[[j]] <- plot_ly(z = z0, x=time, y = u, type = "heatmap", colorscale = list(c(0,'#FFFFFAFF'), c(1,'#FF0000FF')),
                             showscale =FALSE) %>%
            layout(title = y_var, yaxis = list(title = xlab),xaxis = list(title = tlab))
        }else if(type[j]=="3Dsurface"){

          z0 <- Y@B[[j]]%*%Y@C[[j]]
          axx <-axy <-axz <- list(
            gridcolor="rgb(180, 180, 180)",
            zerolinecolor="rgb(255,255,255)"
          )
          axx$title <- ifelse(is.null(tlab),"time",tlab)
          axy$title <- ifelse(is.null(xlab),"x",xlab)
          axz$title <- ifelse(is.null(ylab),paste("Variable",j),ylab[j])
         Pl[[j]] <- plot_ly(z = z0, x=time, y = u, colorscale = list(c(0,'#FFFFFAFF'), c(1,'#FF0000FF'))) %>%
            layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))%>%
            add_surface(showscale=FALSE)

        }else if (type[j]=="3Dline"){

          D0 <- tibble::as_tibble(data.frame(z=c(Y@B[[j]]%*%Y@C[[j]])))
          D0$time <- rep(1:N,each=nrow(Y@grid[[j]]))
          D0$x <- rep(1:nrow(Y@grid[[j]]),ncol(Y@C[[j]]))
          axx <-axy <-axz <- list(
            gridcolor="rgb(180, 180, 180)",
            zerolinecolor="rgb(255,255,255)"
          )
          axx$title <- ifelse(is.null(tlab),"x",tlab)
          axy$title <- ifelse(is.null(xlab),"time",xlab)
          axz$title <- ifelse(is.null(ylab),paste("Variable",var),ylab[j])
          Pl[[j]]<- D0 %>%
            group_by(time) %>%
            plot_ly(x=~time,z=~z,y=~x, type = 'scatter3d', mode = 'lines', color = ~z,
                    line = list(width = 4), colors=c("#FFFFFAFF","#FF0000FF")) %>%
            layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz)) %>% hide_colorbar()


        }

      }else{

        y <- tibble::as_tibble(data.frame(y=c(Y@B[[j]]%*%Y@C[[j]])))
        y$time <- as.factor(rep(time,each=nrow(Y@grid[[j]])))
        y$x_1 <- rep(Y@grid[[j]][,1],ncol(Y@C[[j]]))
        y$x_2 <- rep(Y@grid[[j]][,2],ncol(Y@C[[j]]))
        Pl[[j]] <- ggplotly(ggplot(y,aes(x_2,x_1,fill=y,frame=time))+geom_tile()+scale_fill_distiller(palette = "RdYlBu")+theme_ipsum())

        }


      }

  }else{


    for(j in 1:length(vars)){

      if(ncol(Y@grid[[vars[j]]])==1){

        if(is.na(type[j])==TRUE||type[j]=="line"){

          y <- tibble::as_tibble(data.frame(y=c(Y@B[[vars[j]]]%*%Y@C[[vars[j]]])))
          y$time <- as.factor(rep(time,each=nrow(Y@grid[[vars[j]]])))
          y$x <- rep(1:nrow(Y@grid[[vars[j]]]),ncol(Y@C[[vars[j]]]))
          if(ylab=="y") y_var <- paste("Variable",vars[j]) else y_var <- ylab[vars[j]];
          Pl[[j]] <- y %>%
            group_by(time) %>%
            plot_ly(x=~x,y=~y) %>%
            add_lines(color = ~time,colors=c("lightsteelblue","royalblue4"),showlegend=FALSE) %>%
            layout(title=main,yaxis = list(title = ylab),xaxis = list(title = xlab))


        }else if(type[j]=="heatmap"){

          z0 <- Y@B[[vars[j]]]%*%Y@C[[vars[j]]]
          if(ylab=="y") y_var <- paste("Variable",vars[j]) else y_var <- ylab[vars[j]];
          Pl[[j]] <- plot_ly(z = z0, x=time, y = u, type = "heatmap", colorscale = list(c(0,'#FFFFFAFF'), c(1,'#FF0000FF')),
                             showscale =FALSE) %>%
            layout(title = y_var, yaxis = list(title = xlab),xaxis = list(title = tlab))
        }else if(type[j]=="3Dsurface"){

          z0 <- Y@B[[vars[j]]]%*%Y@C[[vars[j]]]
          axx <-axy <-axz <- list(
            gridcolor="rgb(180, 180, 180)",
            zerolinecolor="rgb(255,255,255)"
          )
          axx$title <- ifelse(is.null(tlab),"time",tlab)
          axy$title <- ifelse(is.null(xlab),"x",xlab)
          axz$title <- ifelse(is.null(ylab),paste("Variable",vars[j]),ylab[vars[j]])
          Pl[[j]] <- plot_ly(z = z0, x=time, y = u, colorscale = list(c(0,'#FFFFFAFF'), c(1,'#FF0000FF'))) %>%
            layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))%>%
            add_surface(showscale=FALSE)

        }else if (type[j]=="3Dline"){

          D0 <- tibble::as_tibble(data.frame(z=c(Y@B[[vars[j]]]%*%Y@C[[vars[j]]])))
          D0$time <- rep(1:N,each=nrow(Y@grid[[vars[j]]]))
          D0$x <- rep(1:nrow(Y@grid[[vars[j]]]),ncol(Y@C[[vars[j]]]))
          axx <-axy <-axz <- list(
            gridcolor="rgb(180, 180, 180)",
            zerolinecolor="rgb(255,255,255)"
          )
          axx$title <- ifelse(is.null(tlab),"x",tlab)
          axy$title <- ifelse(is.null(xlab),"time",xlab)
          axz$title <- ifelse(is.null(ylab),paste("Variable",var),ylab[vars[j]])
          Pl[[j]]<- D0 %>%
            group_by(time) %>%
            plot_ly(x=~time,z=~z,y=~x, type = 'scatter3d', mode = 'lines', color = ~z,
                    line = list(width = 4), colors=c("#FFFFFAFF","#FF0000FF")) %>%
            layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz)) %>% hide_colorbar()


        }

      }else{

        y <- tibble::as_tibble(data.frame(y=c(Y@B[[j]]%*%Y@C[[j]])))
        y$time <- as.factor(rep(time,each=nrow(Y@grid[[j]])))
        y$x_1 <- rep(Y@grid[[j]][,1],ncol(Y@C[[j]]))
        y$x_2 <- rep(Y@grid[[j]][,2],ncol(Y@C[[j]]))
        Pl[[j]] <- ggplotly(ggplot(y,aes(x_2,x_1,fill=y,frame=time))+geom_tile()+scale_fill_distiller(palette = "RdYlBu")+theme_ipsum())

      }


    }



  }
  if("3Dsurface"%in%type || "3Dline"%in%type){

    print(Pl)

  }else if("3Dsurface"%in%type==FALSE && subplot!=TRUE || "3Dline"%in%type==FALSE && subplot!=TRUE){

    print(Pl)

  }else{

    print(subplot(Pl))

  }
}
  #p <- length(Y@coefs)
  #N <- ncol(Y@coefs[[1]])
  #time <- colnames(Y@coefs[[1]])
  #Pl <- list()
  #if(is.null(vars)){
  #  for(i in 1:p) {
  #    if(ncol(Y@grid[[i]])==1) {
  #      if(type1=="line"){
  #        y <- as.tbl(data.frame(y=c(Y@B[[i]]%*%Y@coefs[[i]])))
  #        y$time <- as.factor(rep(time,each=nrow(Y@grid[[i]])))
#           y$x <- rep(rownames(temperature_X),length = nrow(Y@B[[i]])*ncol(Y@coefs[[i]]))
#           if(ylab=="y") y_var <- paste("Variable",i) else y_var <- ylab[i];
#           Pl[[i]] <- y %>%
#           group_by(time) %>%
#           plot_ly(x=~x,y=~y,line=list(width=lwd)) %>%
#           add_lines(color = ~time,colors=c("lightsteelblue","royalblue4"),showlegend=FALSE) %>%
#           layout(title="(A) Intraday Temperature Curves",font=list(size=fmain),margin=list(t=50),yaxis = list(title = "Temperature (C)",gridcolor=toRGB("gray50"),titlefont=list(size=ftitley),tickfont=list(size=fticky)),xaxis = list(title = "Local Standard Time",gridcolor=toRGB("gray50"),titlefont=list(size=ftitley),tickfont=list(size=fticky),ticktext=as.list(temp_times),tickvals=as.list(0:23)))
#           orca(p = Pl[[1]], file = "./mfssa_plots/temp_curves.pdf",width = 600, height = 600)
#         }else if(type1=="heatmap"){
#
#
#           if(is.null(ylab)) y_var <- paste("Variable",i) else y_var <- ylab[i];
#           z0 <- as.matrix(Y@B[[i]]%*%Y@coefs[[i]])
#           Pl[[i]] <- plot_ly(z = ~z0, x= ~time, y = ~unique(Y@grid[[i]][,"x"]), type = "heatmap", colorscale = list(c(0,'#FFFFFAFF'), c(1,'#FF0000FF')),
#                              colorbar=list(title=zlab,lenmode="pixels"))%>%
#             layout(yaxis = list(title = y_var),xaxis = list(title = tlab))
#
#
#
#         }else if(type1=="3Dsurface"){
#
#
#           # to be added later
#
#
#
#         }else{
#
#
#
#           # 3Dline to be added later
#
#
#         }
#       }else{
#         if(type2=="heatmap"){
#           observ<-matrix(data=Y@B[[i]]%*%Y@coefs[[i]][,obs],nrow=length(unique(Y@grid[[i]][,"x"])),ncol=length(unique(Y@grid[[i]][,"y"])))
#           Pl[[i]]=plot_ly(z= ~observ,x= ~unique(Y@grid[[i]][,"x"]), y= ~unique(Y@grid[[i]][,"y"]),type="heatmap",colors=terrain.colors(7),reversescale=T,colorbar=list(title=zlab,lenmode="pixels"))%>%layout(
#           xaxis=list(title=xlab),yaxis=list(title=ylab),title=paste("Observation ",colnames(Y@coefs[[i]])[obs],sep=""))
#         }else{
#             # 3D surface option to view two-dimensional domain observations.
#           observ<-matrix(data=Y@B[[i]]%*%Y@coefs[[i]][,obs],nrow=length(unique(Y@grid[[i]][,"x"])),ncol=length(unique(Y@grid[[i]][,"y"])))
#           Pl[[i]]=plot_ly(z= ~observ,x= ~unique(Y@grid[[i]][,"x"]), y= ~unique(Y@grid[[i]][,"y"]))%>%
#             layout(scene = list(xaxis = list(title=xlab), yaxis = list(title=ylab), zaxis = list(title=zlab)))%>%
#             add_surface(showscale=FALSE)
#
#         }
#
#       }
#
#     }
#   }else{
#
#     for(i in 1:length(vars)){
#       var=vars[i]
#     if(ncol(Y@grid[[var]])==1){
#
#
#
#       if(type1=="line"){
#         y <- as.tbl(data.frame(y=c(Y@B[[i]]%*%Y@coefs[[i]])))
#         y$time <- as.factor(rep(time,each=nrow(Y@grid[[i]])))
#         y$x <- rep(rownames(temperature_X),length = nrow(Y@B[[i]])*ncol(Y@coefs[[i]]))
#         if(ylab=="y") y_var <- paste("Variable",i) else y_var <- ylab[i];
#         Pl[[i]] <- y %>%
#           group_by(time) %>%
#           plot_ly(x=~x,y=~y,line=list(width=lwd)) %>%
#           add_lines(color = ~time,colors=c("lightsteelblue","royalblue4"),showlegend=FALSE) %>%
#           layout(yaxis = list(title = "Celcius",titlefont=list(size=ftitley),tickfont=list(size=fticky)),xaxis = list(title = "Time",titlefont=list(size=ftitley),tickfont=list(size=fticky)))
#
#
#               }else if(type1=="heatmap"){
#
#
#                 if(is.null(ylab)) y_var <- paste("Variable",var) else y_var <- ylab[var];
#                 z0 <- as.matrix(Y@B[[var]]%*%Y@coefs[[var]])
#                 Pl[[i]] <- plot_ly(z = ~z0, x= ~time, y = ~unique(Y@grid[[var]][,"x"]), type = "heatmap", colorscale = Earth,colorbar=list(title=zlab,lenmode="pixels"))%>%
#                   layout(yaxis = list(title = y_var),xaxis = list(title = tlab))
#
#
#
#               }else if(type1=="3Dsurface"){
#
#
#         # to be added later
#
#
#
#               }else{
#
#
#
#         # 3Dline to be added later
#
#
#              }
#
#
#     }else{
#      #### This is the one we want
#       if(type2=="heatmap"){
#               observ<-matrix(data=Y@B[[var]]%*%Y@coefs[[var]][,obs],nrow=length(unique(Y@grid[[var]][,"x"])),ncol=length(unique(Y@grid[[var]][,"y"])))
#               Pl[[i]]=plot_ly(z= ~observ,x= ~unique(Y@grid[[var]][,"x"]), y= ~unique(Y@grid[[var]][,"y"]),type="heatmap",colors=terrain.colors(7),reversescale=T,colorbar=list(title="",lenmode="pixels"))%>%
#                 layout(xaxis=list(title="",ticks="",tickvals=""),yaxis=list(title="",ticks="",tickvals=""),title="NDVI Image Taken September 14, 2009",font=list(size=12))
#
#                 #layout(xaxis=list(side="bottom",title=xlab,titlefont=list(ftitlex)),xaxis2=list(side="top",overlaying="x",tickfont=list(size=ftickx),tickfont=list(ftickx)),yaxis=list(title=ylab,titlefont=list(ftitley),tickfont=list(fticky)),title="NDVI Image Taken September 14, 2009",font=list(size=fmain))
#
#
#               }else{
#         # 3D surface option to view two-dimensional domain observations.
#               observ<-matrix(data=Y@B[[var]]%*%Y@coefs[[var]][,obs],nrow=length(unique(Y@grid[[var]][,"x"])),ncol=length(unique(Y@grid[[var]][,"y"])))
#               Pl[[i]]=plot_ly(z= ~observ,x= ~unique(Y@grid[[var]][,"x"]), y= ~unique(Y@grid[[var]][,"y"]))%>%
#                 layout(scene = list(xaxis = list(title=xlab), yaxis = list(title=ylab), zaxis = list(title=zlab)))%>%
#                 add_surface(showscale=FALSE)
#
#             }
#
#
#     }
#
#
#   }
# }
#
#   Pl2 <- subplot(Pl)
#
#   print(Pl2)



