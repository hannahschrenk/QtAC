#-------------------------------------------------------
#
# QtAC
# Autor: Nico Schreiber
# mod_history:
#
#   DATE                    UPDATE
# 13.06.2019
# 25.06.2019
# 02.07.2019        fixing documentation, modify function names
# 09.07.2019        fixing documentation of functions
# 22.07.2019        version 1.0 of the package
# 08.09.2020	    polishing functions and descriptions
#
#-------------------------------------------------------
# Packages needed:

 # library("pracma")
 # library("rJava")
 # library("igraph")
 # library("rgl")

# "Perl" is necessary as well. You can download it from https://www.perl.org/get.html.

#-------------------------------------------------------
# OVERVIEW ABOUT THE FUNCTIONS:

#                           Input                        		 				Output
# TXT.reader              TXT file                    							data array
# QtAC                    data array                  							list of adjacency and significance matrices
# Signfactor              list of adjacency and significance matrices   				list of adjacency and significance matrices
# writeAdjSgn             list of adjacency and significance matrices   				CSV files containing an adjacency and the corresponding significance matrix respectively
# maturation              list of adjacency and significance matrices   				dataframe containing the three systemic variables of the adjacency matrices
# network                 list of adjacency and significance matrices, number   			network plot of the selected adjacency matrix
# 2dplot                  dataframe containing the systemic variables   			        plot of all or a selected systemic variable w.r.t. time
# 2dmixplot               dataframe containing the systemic variables, names of two variables   	2-dimensional plot of two selected systemic variables w.r.t. each other
# 3dplot		  dataframe containing the systemic variables					3-dimensional plot of the systemic variables w.r.t. each other
#
#-------------------------------------------------------
#' TXT-reader
#'
#' This function is used to import the data in R. The data should be in a tab-seperated file with or without column/row names. Columns should contain time points, rows the system's components.
#' @export
#' @param filename path of the file you want to import
#' @param col_names logical operator. TRUE if the file contains column names, FALSE if it does not
#' @param row_names logical operator. TRUE if the file contains row names, FALSE if it does not
#' @return data array. If no column or row names were given, column or row names of the form "t_ " and "C_ " respectively are added.
# @examples
# data <- QtAC.TXT.reader('example_file.txt', col_names=FALSE, row_names=TRUE)

QtAC.TXT.reader <- function(filename,col_names=FALSE,row_names=FALSE){

Data0 <- read.table(filename, check.names = FALSE, sep = "\t")
d <- dim(Data0)
number_col <- d[2]

if (col_names == FALSE){
 if (row_names == FALSE){
    Data <- read.table(filename, col.names = paste0("t", seq(1,number_col)), check.names = FALSE, sep = "\t")
    num_sp <- seq(1,nrow(Data))
    names_sp <- paste("C",num_sp,sep="_")
    rownames(Data) <- names_sp
 }
 if (row_names == TRUE){
    Data <- read.table(filename, col.names = paste0("t", seq(0,number_col-1)), check.names = FALSE, sep = "\t")
    rownames(Data) <- Data[,1]
    Data <- Data[,-1]
  }
 }

if (col_names == TRUE){

 Data <- read.table(filename, header = TRUE, check.names = FALSE, sep = "\t")

 if (row_names == FALSE){
    num_sp <- seq(1,nrow(Data))
    names_sp <- paste("C",num_sp,sep="_")
    rownames(Data) <- names_sp
    colnames(Data) <- Data[1,]
    Data <- Data[-1,]
 }
 if (row_names == TRUE){
    rownames(Data) <- Data[,1]
    colnames(Data) <- Data[1,]
    Data <- Data[,-1]
    Data <- Data[-1,]
  }
 }
return(Data)
}

#--------------------------------------------------------
#' QtAC (Main function)
#'
#' This function calculates the transfer entropy between two species each for shifting time windows of fixed length.
#' The output is a list of adjacency matrices and the corresponding significance matrices.
#' @export
#' @param Data data array containing time series of the system's components' abundance data
#' @param num_timepoints length of the time windows of abundance data serving as basis of the transfer entropy calculations
#' @param JavaPath path of the file "infodynamics.jar"
#' @param num_PermCheck number of surrogate samples to bootstrap to generate the distribution in the significance test
#' @param k embedding length of destination past history to consider
#' @param k_tau embedding delay for the destination variable
#' @param l embedding length of source past history to consider
#' @param l_tau embedding delay for the source variable
#' @param delay time lag between last element of source and destination next value
#' @return list of two lists (adjacency and corresponding significance matrices)
# @examples
# result_mtx <- QtAC(Data,num_timepoints=6,JavaPath = "D:/Users/max.mustermann/Desktop/infoDynamics.jar",num_PermCheck=500)

QtAC <- function(Data, num_timepoints = 5,JavaPath,num_PermCheck=1000L,k=1L,k_tau=1L,l=1L,l_tau=1L,delay=1L){


  Data1 <- t(Data)


  result_mtx <- list()
  num_timesteps <- dim(Data1)[1]
  num_species <- dim(Data1)[2]

  if(num_timepoints < 5 || num_timepoints > num_timesteps){
    print('num_timepoints has to be between 5 and the length of timepoints in the data')

   } else {
    DataSet <<- .QtAC.split(num_timepoints,num_timesteps,Data1)      #Splitting & Extending

    # Calculation of adjacency and significance matrices:

    adjacent_matrices <- list()
    significance_matrices <- list()


    for (num_dataset in 1:length(DataSet)){
      adj_sign <- .QtAC.Kraskov(DataSet[[num_dataset]],num_species,num_PermCheck,k,k_tau,l,l_tau,delay, rownames(Data),JavaPath = JavaPath)
      adjacent_matrices[[num_dataset]] <- adj_sign[[1]]
      significance_matrices[[num_dataset]] <- adj_sign[[2]]

      Sys.time()     # just for time analysis

      result_mtx[[1]] <- adjacent_matrices
      result_mtx[[2]] <- significance_matrices
      }
  }
  listnames <- colnames(Data)
  names(result_mtx[[1]]) <- listnames[seq(num_timepoints,num_timesteps,1)]
  names(result_mtx[[2]]) <- listnames[seq(num_timepoints,num_timesteps,1)]
  return(result_mtx)
}

#-----------------------------------------------------------
#--------------------------SPLIT------------------------
# This function splits up the data into time windows of predefined length (each of which shifted by one) and extends them by means of a spline if num_timesteps < 15.

.QtAC.split <- function(num_timepoints,num_timesteps,Data){

  DataSet <- list()
  num_DataSet <- 1 + num_timesteps-num_timepoints
  for (k in 1:num_DataSet){
    Submatrix <- Data[k:(k+num_timesteps-num_DataSet),]
    class(Submatrix) <- "double"
    DataSet[[k]] <- Submatrix

  }

  #extend the data if it is to small:

  if(num_timepoints < 15){
    DataSet1 <- list()
    for(row in 1:num_DataSet){
      DataSet1[[row]] <- .QtAC.extend(num_timepoints,DataSet[[row]])
    }
    DataSet <- DataSet1
  }

  return(DataSet)
}

#--------------------------EXTEND------------------------
# This function is used to triple the number of data points with a piecewise cubic spline if num_timesets < 15.

.QtAC.extend <- function(num_timepoints,matrix_step){

  n <- 3*num_timepoints
  tmp <- pracma::zeros(n,ncol(matrix_step))
  ydata <- seq(1,num_timepoints,by=(num_timepoints-1)/(n-1))
  for(col in 1:ncol(matrix_step)){
    tmp[,col] <- pracma::pchip(1:num_timepoints,matrix_step[,col],ydata)
  }
  colnames(tmp) <- colnames(matrix_step)
  return(tmp)
}

#--------------------------Kraskov-------------------------
# This function computes the transfer entropy between two species and stores it in an adjacency matrix. The significance of the results is computed as well and stored in a significance matrix. The return will be a list of two elements, the adjacency and the significance matrix.

.QtAC.Kraskov <- function(Submatrix,num_species,num_PermCheck,k,k_tau,l,l_tau,delay, names_sp,JavaPath){

  rJava::.jinit()
  rJava::.jaddClassPath(JavaPath)
  Submtx_java <- rJava::.jarray(Submatrix, dispatch = TRUE)
  teCalc<-rJava::.jnew("mtinfodynamics/RunTransferEntropyCalculatorKraskov")
  rJava::.jcall(teCalc,"V","initialise",k,k_tau,l,l_tau,delay,num_PermCheck)
  rJava::.jcall(teCalc,"V","runTEKraskov",Submtx_java)
  result_Adj <-  rJava::.jcall(teCalc,"[[D",method = "getResults")
  result_Sgn <- rJava::.jcall(teCalc,"[[D",method = "getSigMatx")


  vector_list_adj <- lapply(result_Adj,function(mat) rJava::.jevalArray(mat, simplify = TRUE))
  vector_list_sgn <- lapply(result_Sgn,function(mat) rJava::.jevalArray(mat, simplify = TRUE))

  Adj <- do.call(cbind, vector_list_adj)
  Sgn <- do.call(cbind, vector_list_sgn)

  Adj[Adj<0] <- 0

  colnames(Adj) <- names_sp
  colnames(Sgn) <- names_sp

  adj_sign <- list(Adj,Sgn)

  return(adj_sign)

}


#-----------------------------------------------------------
#' Signfactor
#'
#' This function sets all entries in the adjacency matrices to 0 whose p-value is above the predefined significance level.
#' @export
#' @param result_mtx list of adjacency and significance matrices
#' @param signfac significance level
#' @return list of two lists (adjacency matrix containing only significant transfers and significance matrices)
# @examples
# result_mtx_new <- QtAC.Signfactor(result_mtx,Signfac = 0.05)

QtAC.Signfactor <- function(result_mtx,signfac = 0.1){

  adjacent_matrices <- result_mtx[[1]]
  significance_matrices <- result_mtx[[2]]
  transformed_mtx <- list()
  adj_mtx <- list()

  for (num_mtx in 1:length(adjacent_matrices)){

    sign_matrix <- significance_matrices[[num_mtx]]
    adj_matrix <- adjacent_matrices[[num_mtx]]

    adj_matrix[sign_matrix > signfac] = 0

    adj_mtx[[num_mtx]] <- adj_matrix
  }
  transformed_mtx[[1]] <- adj_mtx
  transformed_mtx[[2]] <- significance_matrices

  names(transformed_mtx[[1]]) <- names(result_mtx[[1]])
  names(transformed_mtx[[1]]) <- names(result_mtx[[2]])

  return(transformed_mtx)
}


#----------------------------------------------------------
#' WriteAdjSgn
#'
#' This function writes a selected adjacency and significance matrix to CSV files.
#' @export
#' @param result_mtx list of adjacency and significance matrices
#' @param num_mtx number of the adjacency and significance matrix you want to write to a CSV file
#' @param file_name_adj name of the file to which the adjacency matrix is written
#' @param file_name_sgn name of the file to which the significance matrix is written
#' @return two CSV files, one containing the selected adjacency, one the significance matrix
# @examples
# QtAC.WriteAdjSgn(result_mtx,2,file_name_adj = "adj_mtx",file_name_sgn = "sgn_mtx")

QtAC.WriteAdjSgn <- function(result_mtx,num_mtx,file_name_adj = "adjacency_matrix",file_name_sgn = "significance_matrix"){

  adjacent_matrices <- result_mtx[[1]]
  significance_matrices <- result_mtx[[2]]


  write.table(as.matrix(adjacent_matrices[[num_mtx]]), file=paste(file_name_adj,"_",num_mtx,".csv", sep=""), sep = "\t", row.names=F,quote=F)
  write.table(as.matrix(significance_matrices[[num_mtx]]), file=paste(file_name_sgn,num_mtx,".csv", sep=""), sep = "\t", row.names=F,quote=F)

}


#---------------------------------------------------------------
#---------------------------Potential---------------------------
# This function calculates the potential of a network.

.QtAC_potential <- function(mtx){

  mtx <- round(mtx,digits = 3)

  if(all(mtx==0)){
    potential <- 0

  }
  else{
    length <- ncol(mtx)

    norm_mtx <- mtx/sum(mtx)
    b = 0
    for(source in 1:length){
      a = 0
      for(dest in 1:length){
        if(norm_mtx[source,dest] != 0){
          a <- a + norm_mtx[source,dest] * log2(norm_mtx[source,dest])
        }
      }
      b <- b+a
    }
    potential <- -b * sum(mtx)

  }

  return(potential)

}
#-------------------------Connectedness-------------------------
# This function calculates the connectedness of a network.

.QtAC.connectedness <- function(mtx){

  mtx <- round(mtx,digits = 3)

  length <- ncol(mtx)
  b <- 0

  for(i in 1:length){
    a <-0
    for(j in 1:length){

      colval <- sum(mtx[,j])
      rowval <- sum(mtx[i,])

      if(mtx[i,j] != 0){
        a <- a + mtx[i,j] * log2(mtx[i,j]*sum(mtx)/(colval*rowval))
      }
    }
    b <- b+a
  }

  connectedness <- b
  return(connectedness)

}

#---------------------------Resilience-----------------------
# This function calculates the resilience of a network.


.QtAC.resilience <- function(mtx, norm = res_norm, stand = res_stand){

  mtx <- round(mtx,digits = 3)

  length <- ncol(mtx)
  L_out <- pracma::zeros(length,length)
  L_in <- pracma::zeros(length,length)

if(norm == "continous"){

  for(i in 1:length){
   for(j in 1:length){
     if(i==j){

       L_out[i,j] <- sqrt(sum(mtx[i,]))
     }
     else if(sum(mtx[i,]) != 0){

       L_out[i,j] <- -mtx[i,j]/sqrt(sum(mtx[i,]))

     }
   }
  }



  for(i in 1:length){
    for(j in 1:length){
      if(i==j){

        L_in[i,j] <- sqrt(sum(mtx[,j]))
      }
      else if(sum(mtx[,j]) != 0){

        L_in[i,j] <- -mtx[i,j]/sqrt(sum(mtx[,j]))

      }
    }
  }

}

if(norm == "symmetric"){

  for(i in 1:length){
   for(j in 1:length){
     if(i==j && sum(mtx[i,])!=0){

       L_out[i,j] <- 1 - mtx[i,j]/sum(mtx[i,])
     }
     if(i != j && sum(mtx[i,]) != 0 && sum(mtx[j,]) != 0){

       L_out[i,j] <- -mtx[i,j]/sqrt(sum(mtx[i,])*sum(mtx[j,]))

     }
   }
  }



  for(i in 1:length){
    for(j in 1:length){
      if(i==j && sum(mtx[,i])!=0){

        L_in[i,j] <- 1 - mtx[i,j]/sum(mtx[,i])
      }
      if(i != j && sum(mtx[,j]) != 0 && sum(mtx[,i]) != 0){

        L_in[i,j] <- -mtx[i,j]/sqrt(sum(mtx[,j])*sum(mtx[,i]))

      }
    }
  }

}



if(stand == "maxweight"){


  if(!(all(mtx==0))){

    b <- max(mtx)

    L_out <- (1/sqrt(b)) * L_out
    L_in <- (1/sqrt(b)) * L_in

    }

}

if(stand == "maxweight2"){


  if(!(all(mtx==0))){

    b <- max(mtx)

    L_out <- (1/b) * L_out
    L_in <- (1/b) * L_in

    }

}

if(stand == "nodes"){


  if(!(all(mtx==0))){

    b <- max(mtx)

    L_out <- (sqrt(length-1)/length) * L_out
    L_in <- (sqrt(length-1)/length) * L_in

    }

}


if(stand == "maxweightnodes"){


  if(!(all(mtx==0))){

    b <- max(mtx)

    L_out <- (sqrt(length-1)/length) * 1/sqrt(b) * L_out
    L_in <- (sqrt(length-1)/length) * 1/sqrt(b) * L_in

    }

}

  D_out <- round(abs(Re(eigen(L_out,only.values = TRUE)[[1]])),10)
  D_in <- round(abs(Re(eigen(L_in,only.values = TRUE)[[1]])),10)

  if(sum(D_out)!=0){

   d_out <- min(D_out[D_out != 0])

  } else {
    d_out <- 0
  }
  if(sum(D_in)!=0){

    d_in <- min(D_in[D_in != 0])

  } else {
    d_in <- 0
  }
  resilience <- min(d_in,d_out)

  return(resilience)

}
#------------------------------------------------------------------
#' maturation
#'
#' This function computes the three systemic variables (potential, connectedness, and resilience) of each adjacency matrix.
#' @export
#' @param result_mtx list of adjacency matrices and significance matrices
#' @param res_norm normalization variant of the Laplacian matrices ("continous", "symmetric"). If res_norm = "symmetric", \eqn{L_{out} = c \cdot D^{-\frac{1}{2}}_{out} \left(D_{out}- A \right)} and \eqn{L_{in} = c \cdot \left(D_{in} - A \right) D^{-\frac{1}{2}}_{in}}. If res_norm = "symmetric", \eqn{L_{out} = c \cdot D^{-\frac{1}{2}}_{out} (D_{out}-A) D^{-\frac{1}{2}}_{out}} and \eqn{L_{in} =  c\cdot D^{-\frac{1}{2}}_{in} (D_{in}-A)D^{-\frac{1}{2}}_{in}}.
#' @param res_stand standardization constant c of the Laplacian matrices ("none", "maxweight", "maxweight2", "nodes", "maxweightnodes"). Let \eqn{N} be the number of nodes of the underlying graph and \eqn{M} its maximal edge weight. If res_stand = "none", \eqn{c = 1}. If res_stand = "maxweight", \eqn{c = \frac{1}{\sqrt{M}}}. If res_stand = "maxweight2", \eqn{c = \frac{1}{M}}. If res_stand = "nodes", \eqn{c = \frac{\sqrt{N-1}}{N}}. If res_stand = "maxweightnodes", \eqn{c = \frac{\sqrt{N-1}}{N \cdot \sqrt{M}}}.
#' @return dataframe containing the three systemic variables (potential, connectedness, and resilience) of each adjacency matrix
# @examples
# Mat <- QtAC.maturation(result_mtx, res_norm = "continous", res_stand = "maxweightnodes")

QtAC.maturation <- function(result_mtx, res_norm = "continous", res_stand = "maxweight"){

  Adjacent_matrices <- result_mtx[[1]]
  length <- length(Adjacent_matrices)

  potential <- pracma::zeros(length,1)
  connectedness <- pracma::zeros(length,1)
  resilience <- pracma::zeros(length,1)

  for(num in 1:length){

   potential[num] <- .QtAC_potential(Adjacent_matrices[[num]])
   connectedness[num] <- .QtAC.connectedness(Adjacent_matrices[[num]])
   resilience[num] <- .QtAC.resilience(Adjacent_matrices[[num]], norm = res_norm, stand = res_stand)

  }

  maturation <- data.frame(potential,connectedness,resilience)
  rownames(maturation) <- names(result_mtx[[1]])
  return(maturation)

}

#----------------------------------------------------------------
#' network
#'
#' This function plots a selected adjacency matrix as network.
#' @export
#' @param result_mtx list of adjacency matrices and significance matrices
#' @param num_mtx number of the adjacency matrix you want to plot
#' @param edge_label If edge_label = TRUE, the weight of the edges are plotted next to the edges.
#' @param dec number of decimal digits in the edge labels
#' @param layout layout format ("circle","star","fruchterman.reingold","grid","nicely")
#' @param edge_width muliplicator for the width of the edges
#' @param arrow_width muliplicator for the width of the arrows
#' @param col_node color of the vertices
#' @param col_edge color of the edges
#' @param vertex_label If vertex_label = "short", the first 3 letters of the components' names will be used as vertex labels, if vertex_label = "long", the whole names will be used as vertex labels. Via vertex_label =  c(...), customized names can be used as vertex labels.
#' @param Save If Save = TRUE, the network will be saved as a PNG file.
#' @param filename If Save = TRUE, the network will be saved in a file called filename.
#' @return network plot and, if Save = TRUE, a PNG file containing the plot
# @examples
# QtAC.network(result_mtx,3,edge_label = TRUE, dec = 3, layout = "nicely", edge_width = 30, arrow_width = 30,col_node = "red", vertex_label = "short", Save = FALSE)

QtAC.network <- function(result_mtx,num_mtx,edge_label = FALSE,dec = 2, layout = "circle",edge_width = 3,arrow_width = 5,col_node = "palegreen3",col_edge = "steelblue3",vertex_label = "short", Save = FALSE,filename = "network"){

  adj_mtx <- result_mtx[[1]][[num_mtx]]

  graph_adj_mtx <- igraph::graph.adjacency(adj_mtx, mode="directed", weighted=TRUE)

  times <- names(result_mtx[[1]])

  title <- paste("Network of Adjacency matrix ",times[num_mtx])


   if(layout=="circle"){


     layout_graph <- igraph::layout_in_circle(graph_adj_mtx)

   }

   else if(layout=="star"){

     layout_graph <- igraph::layout_as_star(graph_adj_mtx)

   }

   else if(layout == "fruchterman.reingold"){

     layout_graph <- igraph::layout.fruchterman.reingold(graph_adj_mtx)

   }

   else if(layout == "grid"){

     layout_graph <- igraph::layout_on_grid(graph_adj_mtx)

    }

  else if(layout == "nicely"){

    layout_graph <- igraph::layout_nicely(graph_adj_mtx)

  }

   else {

     return("Layout not available")
   }

  #--------------------------------------------------

  if(edge_label ==TRUE){
  edge_label <- round(igraph::E(graph_adj_mtx)$weight,dec)
  }
  else edge_label = NULL

  #----------------------------------------------------

  if(vertex_label == "long"){

    vertex_label <- colnames(adj_mtx)
  }
  else if(vertex_label == "short"){

    names_long <- colnames(adj_mtx)
    vertex_label <- character(length(names_long))

    for (name in 1:length(names_long)){

    vertex_label[name] <- paste(strsplit(names_long[name],"")[[1]][1],strsplit(names_long[name],"")[[1]][2],strsplit(names_long[name],"")[[1]][3],sep = "")

    }
  }

  #----------------------------------------------------

   if(Save){
    filename <- paste(paste(filename,num_mtx,sep = '_'),'.png',sep='')
    png(filename,width = 800 ,height= 800,units = "px",type = "cairo")
    igraph::plot.igraph(graph_adj_mtx,layout=layout_graph, vertex.label = vertex_label, vertex.label.color = "black",vertex.color = col_node,main= title,edge.label = edge_label, edge.label.font= 2, edge.arrow.size =igraph::E(graph_adj_mtx)$weight/max(igraph::E(graph_adj_mtx)$weight)*arrow_width, edge.width=igraph::E(graph_adj_mtx)$weight/max(igraph::E(graph_adj_mtx)$weight)*edge_width, vertex.label.cex = 1, margin = 0, vertex.label.font=2, edge.curved= 0.1 , edge.color= col_edge, asp= 1)
    dev.off()
   }

  #plot the network, no matter whether it should be saved or not
  igraph::plot.igraph(graph_adj_mtx,layout=layout_graph, vertex.label = vertex_label, vertex.label.color = "black",vertex.color = col_node,main= title,edge.label = edge_label, edge.label.font= 2, edge.arrow.size =igraph::E(graph_adj_mtx)$weight/max(igraph::E(graph_adj_mtx)$weight)*arrow_width, edge.width=igraph::E(graph_adj_mtx)$weight/max(igraph::E(graph_adj_mtx)$weight)*edge_width, vertex.label.cex = 1, margin = 0, vertex.label.font=2, edge.curved= 0.1 , edge.color= col_edge, asp= 1)

  }

#----------------------------------------------------------------
#' 2dplot
#'
#' This function plots graphs of potential, connectedness and resilience with respect to time.
#' @export
#' @param Mat dataframe containing the systemic variables
#' @param prop If prop = NULL, the three systemic variables are plotted w.r.t time in one plots each. If prop = "potential", "connectedness", or "resilience", only the selected systemic variable is plotted w.r.t time.
#' @param time_int vector containing start time, end time, and step size to define the xaxis
#' @param Save If Save = TRUE, the 2D plot will be saved in a PNG file.
#' @param filename If Save = TRUE, the network will be saved in a file called filename.
#' @return 2D plot and, if Save = TRUE, a PNG file containing the plots.
# @examples
# QtAC.2dplot(maturation,prop = "potential", time_int = c(1990,2018,2))

QtAC.2dplot <- function(Mat,prop = NULL, time_int = NULL, Save = FALSE, filename = "2dplot"){

  if(is.null(prop)){
    if(Save){
      filename <- paste(filename,'.png',sep='')
      png(filename,width = 800 ,height= 800,units = "px",type = "cairo")
      #png(filename,units = "px",type = "cairo")
      par(mfrow=c(3,1))
      .QtAC.potentialplot(Mat[,1],time_int,rownames(Mat),FALSE,"")
      .QtAC.connectednessplot(Mat[,2],time_int,rownames(Mat),FALSE,"")
      .QtAC.resilienceplot(Mat[,3],time_int,rownames(Mat),FALSE,"")
      dev.off()

      par(mfrow=c(3,1))
      .QtAC.potentialplot(Mat[,1],time_int,rownames(Mat),FALSE,"")
      .QtAC.connectednessplot(Mat[,2],time_int,rownames(Mat),FALSE,"")
      .QtAC.resilienceplot(Mat[,3],time_int,rownames(Mat),FALSE,"")
    }
    else {
    .QtAC.potentialplot(Mat[,1],time_int,,rownames(Mat), FALSE,"")
    .QtAC.connectednessplot(Mat[,2],time_int,rownames(Mat), FALSE,"")
    .QtAC.resilienceplot(Mat[,3],time_int,rownames(Mat), FALSE,"")
    }
  }

  else if(prop=="potential"){

    vec <- Mat[,1]
    par(mfrow = c(1,1))
    .QtAC.potentialplot(vec,time_int,rownames(Mat),Save,filename)

  }

  else if(prop=="connectedness"){

    vec <- Mat[,2]
    par(mfrow = c(1,1))
    .QtAC.connectednessplot(vec,time_int,rownames(Mat),Save,filename)
  }

  else if(prop=="resilience"){

    vec <- Mat[,3]
    par(mfrow = c(1,1))
    .QtAC.resilienceplot(vec,time_int,rownames(Mat),Save,filename)

  }

  else return("No available property")



}

#------------------------------Potentialplot------------------------------
# This function plots the potential w.r.t time. It is used in QtAC.2dplot.

.QtAC.potentialplot <- function(vec, time_int = NULL,names = seq(1,length(vec),1),Save,filename){


  xaxis_mat <- 0:(length(vec)-1)
  xaxis_mat_spl <- seq(0,length(vec)-1,by = 0.01)


  if(is.null(time_int)){

    xaxis_plot <- 1:(length(vec))
    xaxis_plot_spl <- seq(1,length(vec),by = 0.01)
    start_val <- 1
    end_val <- length(vec)
    timestep <- 1
  }

  else{

    num_timepoints <- (time_int[2]-time_int[1])/time_int[3] + 2 - length(vec)
    start_val <- time_int[1] + time_int[3]*(num_timepoints-1)
    end_val <- time_int[2]

    xaxis_plot <- seq(start_val,end_val,by = time_int[3])
    xaxis_plot_spl <- seq(start_val,end_val,by = ((end_val-start_val)*0.01)/(length(vec)-1))
    timestep <- time_int[3]

  }

  #-------------#
  # splining of the maturation curve:

  vec_ext <- pracma::pchip(xaxis_mat,vec,xaxis_mat_spl)


  xaxis_values <- seq(start_val,end_val, by = timestep)

  # if Save = TRUE the plot will be saved in a png file:

  if(Save){
    filename <- paste(filename,'.png',sep='')
    png(filename,width = 800 ,height= 800,units = "px",type = "cairo")

  plot(xaxis_plot_spl,vec_ext,type = "l",col = "darkgreen",main = paste("Potential"),col.main = "darkblue",xlab = "", ylab = "", lwd = 3, axes = FALSE, frame.plot = TRUE)
  axis(side = 1, at = xaxis_values, labels = names, col.axis = "darkblue")
  axis(side = 2, col.axis = "darkblue")
  mtext("time",side = 1,line = 3, col = "darkblue", font = 2)
  mtext("Potential",side = 2,line = 3, col = "darkgreen", font = 2)
  dev.off()
  }

  plot(xaxis_plot_spl,vec_ext,type = "l",col = "darkgreen",main = paste("Potential"),col.main = "darkblue",xlab = "", ylab = "", lwd = 3, axes = FALSE, frame.plot = TRUE)
  axis(side = 1, at = xaxis_values,labels = names, col.axis = "darkblue")
  axis(side = 2, col.axis = "darkblue")
  mtext("time",side = 1,line = 3, col = "darkblue", font = 2)
  mtext("Potential",side = 2,line = 3, col = "darkgreen", font = 2)

}

#--------------------------Connectednessplot-------------------------
# This function plots the connectedness w.r.t time. It is used in QtAC.2dplot.

.QtAC.connectednessplot <- function(vec, time_int = NULL,names = seq(1,length(vec),1),Save,filename){


  xaxis_mat <- 1:(length(vec))
  xaxis_mat_spl <- seq(1,length(vec),by = 0.01)


  if(is.null(time_int)){

    xaxis_plot <- 1:(length(vec))
    xaxis_plot_spl <- seq(1,length(vec),by = 0.01)
    start_val <- 1
    end_val <- length(vec)
    timestep <- 1

  }

  else{

    num_timepoints <- (time_int[2]-time_int[1])/time_int[3] + 2 - length(vec)
    start_val <- time_int[1] + time_int[3]*(num_timepoints-1)
    end_val <- time_int[2]

    xaxis_plot <- seq(start_val,end_val,by = time_int[3])
    xaxis_plot_spl <- seq(start_val,end_val,by = ((end_val-start_val)*0.01)/(length(vec)-1))
    timestep <- time_int[3]

  }

  #-------------#
  # splining of the maturation curve:

  vec_ext <- pracma::pchip(xaxis_mat,vec,xaxis_mat_spl)


  xaxis_values <- seq(start_val,end_val, by = timestep)


  if(Save){
    filename <- paste(filename,'.png',sep='')
    png(filename,width = 800 ,height= 800,units = "px",type = "cairo")

  plot(xaxis_plot_spl,vec_ext,type = "l",col = "blue",main = paste("Connectedness"),col.main = "darkblue",xlab = "", ylab = "", lwd = 3, axes = FALSE, frame.plot = TRUE)
  axis(side = 1, at = xaxis_values,labels = names, col.axis = "darkblue")
  axis(side = 2, col.axis = "darkblue")
  mtext("time",side = 1,line = 3, col = "darkblue", font = 2)
  mtext("Connectedness",side = 2,line = 3, col = "blue", font = 2)
  dev.off()
  }

  plot(xaxis_plot_spl,vec_ext,type = "l",col = "blue",main = paste("Connectedness"),col.main = "darkblue",xlab = "", ylab = "", lwd = 3, axes = FALSE, frame.plot = TRUE)
  axis(side = 1, at = xaxis_values,labels = names, col.axis = "darkblue")
  axis(side = 2, col.axis = "darkblue")
  mtext("time",side = 1,line = 3, col = "darkblue", font = 2)
  mtext("Connectedness",side = 2,line = 3, col = "blue", font = 2)

}

#-----------------------------Resilienceplot-----------------------------
# This function plots the resilience w.r.t time. It is used in QtAC.2dplot.

.QtAC.resilienceplot <- function(vec, time_int = NULL,names = seq(1,length(vec),1),Save,filename){


  xaxis_mat <- 1:(length(vec))
  xaxis_mat_spl <- seq(1,length(vec),by = 0.01)


  if(is.null(time_int)){

    xaxis_plot <- 1:(length(vec))
    xaxis_plot_spl <- seq(1,length(vec),by = 0.01)
    start_val <- 1
    end_val <- length(vec)
    timestep <- 1

  }

  else{

    num_timepoints <- (time_int[2]-time_int[1])/time_int[3] + 2 - length(vec)
    start_val <- time_int[1] + time_int[3]*(num_timepoints-1)
    end_val <- time_int[2]

    xaxis_plot <- seq(start_val,end_val,by = time_int[3])
    xaxis_plot_spl <- seq(start_val,end_val,by = ((end_val-start_val)*0.01)/(length(vec)-1))
    timestep <- time_int[3]

  }

  #-------------#
  # splining of the maturation curve:

  vec_ext <- pracma::pchip(xaxis_mat,vec,xaxis_mat_spl)


  xaxis_values <- seq(start_val,end_val, by = timestep)

  if(Save){
    filename <- paste(filename,'.png',sep='')
    png(filename,width = 800 ,height= 800,units = "px",type = "cairo")

  plot(xaxis_plot_spl,vec_ext,type = "l",col = "red",main = paste("Resilience"),col.main = "darkblue",xlab = "", ylab = "", lwd = 3, axes = FALSE, frame.plot = TRUE)
  axis(side = 1, at = xaxis_values,labels = names, col.axis = "darkblue")
  axis(side = 2, col.axis = "darkblue")
  mtext("time",side = 1,line = 3, col = "darkblue", font = 2)
  mtext("Resilience",side = 2,line = 3, col = "red", font = 2)
  dev.off()
  }

  plot(xaxis_plot_spl,vec_ext,type = "l",col = "red",main = paste("Resilience"),col.main = "darkblue",xlab = "", ylab = "", lwd = 3, axes = FALSE, frame.plot = TRUE)
  axis(side = 1, at = xaxis_values,labels = names, col.axis = "darkblue")
  axis(side = 2, col.axis = "darkblue")
  mtext("time",side = 1,line = 3, col = "darkblue", font = 2)
  mtext("Resilience",side = 2,line = 3, col = "red", font = 2)

}

#-------------------------------------------------------------------
#' 2dmixplot
#'
#' This function plots two selected systemic variables w.r.t. each other.
#' @export
#' @param Mat data frame containing the systemic variables
#' @param prop1 variable on x-axis ("potential","connectedness","resilience")
#' @param prop2 variable on y-axis ("potential","connectedness","resilience")
#' @param Save If Save = TRUE, the 2D plot will be saved in a PNG file.
#' @param filename If Save = TRUE, the network will be saved in a file called filename.
#' @return 2D plot and, if Save = TRUE, a PNG file containing the plot
# @examples
# QtAC_2dmixplot(Mat,"connectedness","resilience",Save = TRUE, filename = "MIX_plot")

QtAC.2dmixplot <- function(Mat,prop1, prop2,Save = FALSE,filename = paste("2dmixplot", prop1, prop2, sep = "_")){

  if((prop1=="potential")&(prop2=="connectedness")){

    prop_num1 <- 1
    prop_num2 <- 2
    color <- "maroon"

  }

  else if((prop1=="potential")&(prop2=="resilience")){

    prop_num1 <- 1
    prop_num2 <- 3
    color <- "#008080"	 #color: Teal

  }

  else if((prop1=="connectedness")&(prop2=="resilience")){

    prop_num1 <- 2
    prop_num2 <- 3
    color <- "navy"

  }

  else if((prop2=="potential")&(prop1=="connectedness")){

    prop_num1 <- 2
    prop_num2 <- 1
    color <- "maroon"

  }

  else if((prop2=="potential")&(prop1=="resilience")){

    prop_num1 <- 3
    prop_num2 <- 1
    color <- "#008080"	 #color: Teal

  }

  else if((prop2=="connectedness")&(prop1=="resilience")){

    prop_num1 <- 3
    prop_num2 <- 2
    color <- "navy"

  }

  else{
    print(prop1)
    print(prop2)
    return("Property 1 is not avaiable")
  }

  vec1 <- Mat[,prop_num1]
  vec2 <- Mat[,prop_num2]

  vec1_ext <- pracma::pchip(0:(length(vec1)-1),vec1,seq(0,length(vec1)-1,by = 0.1))
  vec2_ext <- pracma::pchip(0:(length(vec2)-1),vec2,seq(0,length(vec2)-1,by = 0.1))

  filename <- paste(filename,".png",sep = "")

  par(mfrow = c(1,1))

  if(Save){
    png(filename,width = 800 ,height= 800,units = "px",type = "cairo")

    plot(vec1_ext,vec2_ext,type = "l",col = color,main = paste(prop1," vs. ", prop2),xlab = prop1, ylab = prop2, lwd = 3)
    text(vec1,vec2,rownames(Mat))
    dev.off()
  }

  plot(vec1_ext,vec2_ext,type = "l",col = color,main = paste(prop1," vs. ", prop2),xlab = prop1, ylab = prop2, lwd = 3)
  text(vec1,vec2,rownames(Mat))
}


#-------------------------------------------------------------
#' 3dplot
#'
#' 3D plot of the three systemic variables w.r.t each other.
#' @export
#' @param Mat data frame containing the three systemic variables
#' @param Mat_points If Mat_points = TRUE, the maturation points are visible.
#' @return 3D plot
# @examples
# QtAC.3dplot(Mat = maturation, Mat_points = TRUE)

QtAC.3dplot <- function(Mat,Mat_points=FALSE){

  steps <- seq(1,dim(Mat)[1],0.01)

  x <- Mat[,1]
  y <- Mat[,2]
  z <- Mat[,3]

  x <- smooth.spline( x, spar=0 )
  y <- smooth.spline( y, spar=0 )
  z <- smooth.spline( z, spar=0 )

  x <- predict(x,steps)$y
  y <- predict(y,steps)$y
  z <- predict(z,steps)$y

  if(Mat_points==TRUE){

    pointsize <- 0.5
  }
  else pointsize <- 0

  rgl::plot3d(Mat[,1], Mat[,2], Mat[,3], pch=19, cex=0.25, size=pointsize,type = "s", col="darkblue",xlab = "Potential",ylab = "Connectedness",zlab = "Resilience")
  rgl::text3d(Mat[,1],Mat[,2], Mat[,3], rownames(Mat))

  rgl::lines3d(x,y,z,col = "red",lwd=6)

}

#-------------------------------------------------------------
