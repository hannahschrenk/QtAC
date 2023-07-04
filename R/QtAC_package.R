#-------------------------------------------------------
#
# QtAC
# Autor: Hannah Schrenk and Nico Schreiber and Carlos Garcia-Perez
# mod_history:
#
#   DATE                    UPDATE
# 13.06.2019
# 25.06.2019
# 02.07.2019        fixing documentation, modify function names
# 09.07.2019        fixing documentation of functions
# 22.07.2019        version 1.0 of the package
# 08.09.2020	      polishing functions and descriptions
# 27.05.2021        adding parameter noise_level, changing descriptions
# 29.11.2021	      changing parameter names
# 07.12.2021 	      version 1.1 of the package
# 04.07.2023        adding functions for estimating AIS, CTE, and local versions of AIS and (C)TE
#                   version 1.2 of the package
#-------------------------------------------------------
# Packages needed:

# library("pracma")
# library("rJava")
# library("igraph")
# library("rgl")

#-------------------------------------------------------
# OVERVIEW ABOUT THE FUNCTIONS:

#                           Input                        		 				                              Output
# TXT.reader              TXT file                    							                              data array
# QtAC                    data array                  							                              list of adjacency and significance matrices
# signfactor              list of adjacency and significance matrices   				                  list of adjacency and significance matrices
# QtAC_CTE                data array                                                              list of adjacency and significance matrices
# signfactor_CTE          list of adjacency and significance matrices                             list of adjacency and significance matrices
# QtAC_AIS                data_array                                                              list of adjacency and significance matrices
# maturation              list of adjacency and significance matrices   				                  data frame containing the three systemic variables of the adjacency matrices
# network                 list of adjacency and significance matrices,                            network plot of the selected adjacency matrix
# 2dplot                  dataframe containing the systemic variables   			                    plot of all or a selected systemic variable w.r.t. time
# 2dmixplot               dataframe containing the systemic variables,                           	2-dimensional plot of two selected systemic variables w.r.t. each other
# 3dplot		              dataframe containing the systemic variables					                    3-dimensional plot of the systemic variables w.r.t. each other
#
#-------------------------------------------------------
#' TXT-reader
#'
#' This function is used to import the data in R. The data should be in a tab-separated file with or without column/row names. Columns should contain time points, rows the system's components.
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
    }
    if (row_names == TRUE){
      rownames(Data) <- Data[,1]
      Data <- Data[,-1]
    }
  }
  return(Data)
}

#-------------------------CTE-------------------------------
#' QtAC_CTE
#'
#' This function calculates the collective transfer entropy of each species for shifting time windows of fixed length or local values.
#' The output is a list of adjacency matrices and the corresponding significance matrices.
#' @export
#' @param data data array containing time series of the system's components' abundance data
#' @param num_timepoints length of the time windows of abundance data serving as basis of the collective transfer entropy estimations
#' @param javapath path of the file "MTinfodynamics.jar"
#' @param noise_level amount of random Gaussian noise added in the estimation
#' @param num_permcheck number of surrogate samples to bootstrap to generate the distribution in the significance test
#' @param k embedding length of destination past history to consider
#' @param k_tau embedding delay for the destination variable
#' @param l embedding length of source past history to consider
#' @param l_tau embedding delay for the source variable
#' @param delay time lag between last element of source and destination next value
#' @param condEmbedDims array of embedding lengths for each conditional variable
#' @param cond_taus array of embedding delays for the conditional variables
#' @param cond_delays array of time lags between last element of each conditional variable and destination next value
#' @param mode Transfer entropy is either averaged over time windows ("average") or local values are estimated ("local"). In the latter case, the value of num_timepoints is ignored.
#' @param save If save=TRUE, the output is saved in a file called filename.
#' @return list of two lists (adjacency and corresponding significance matrices)
# @examples
# result_mtx <- QtAC(Data,num_timepoints=6,JavaPath = "D:/Users/max.mustermann/Desktop/infoDynamics.jar",num_PermCheck=500)

QtAC_CTE <- function(data,num_timepoints = 5,javapath,noise_level = "1e-20",num_permcheck=1000L,k=1L,k_tau=1L,l=1L,l_tau=1L,delay=1L, condEmbedDims = 1L, cond_taus = 1L, cond_delays = 1L, mode = c("average", "local"), save = FALSE, filename = "result_QtAC"){

  Data1 <- t(data)

  result_mtx <- list()
  num_timesteps <- dim(Data1)[1]
  num_species <- dim(Data1)[2]
  names_sp <- colnames(Data1)
  
  result_mtx[[1]] <- list()
  result_mtx[[2]] <- list()

  if(num_timepoints < 5 || num_timepoints > num_timesteps){
    print('num_timepoints has to be between 5 and the length of timepoints in the data')

  	} else if(mode == "local"){ 

	for(t in 1:num_timesteps){
		result_mtx[[1]][[t]] <- matrix(NaN, num_species, num_species-1)
		result_mtx[[2]][[t]] <- matrix(NaN, num_species, num_species-1)
		}

	DataSet <- .QtAC.split(num_timesteps,num_timesteps,Data1)      

        rJava::.jinit()
        rJava::.jaddClassPath(javapath)

	for(i in 1:num_species){

    Submatrix <- DataSet[[1]]	
    destination <- Submatrix[,i]
		Submatrix_temp <- Submatrix[,-i]

		for(j in 1:(num_species-1)){

		  source <- Submatrix_temp[,j]
      mtx_temp <- cbind(source, destination)
      mtx_temp <- rJava::.jarray(mtx_temp, dispatch = TRUE)
		  
			if(j==1){

			  teCalc <- rJava::.jnew("mtinfodynamics/RunTransferEntropyCalculatorKraskov")
  			rJava::.jcall(teCalc,"V","setProperty", "NOISE_LEVEL_TO_ADD", noise_level)
  			rJava::.jcall(teCalc,"V","initialise",k,k_tau,l,l_tau,delay,num_permcheck)
  			rJava::.jcall(teCalc,"V","runTEKraskov",mtx_temp, "local")
			
  			for(t in 0:(num_timesteps-1)){
  			  mtx <- result_mtx[[1]][[t+1]]
  			  rmtx <- rJava::.jcall(teCalc,"[[D",method = "getResultsMtxTimeStep", t)
  			  vector_list_rmtx <- lapply(rmtx, function(mat) rJava::.jevalArray(mat, simplify = TRUE))
  			  rmtx <- do.call(cbind, vector_list_rmtx)
  			  mtx[i,j] <- rmtx[1,2]
  			  result_mtx[[1]][[t+1]] <- mtx
  			  mtx <- result_mtx[[2]][[t+1]]
  			  Sgn <- rJava::.jcall(teCalc,"[[D",method = "getSigMatx")
  			  vector_list_sgn <- lapply(Sgn,function(mat) rJava::.jevalArray(mat, simplify = TRUE))
  			  smtx <- do.call(rbind, vector_list_sgn)
  			  mtx[i,j] <- smtx[1,2]
  			  result_mtx[[2]][[t+1]] <- mtx
  			  }

	   		}

			if(j>1){

			mtx <- cbind(source,destination)
			mtx_java <- rJava::.jarray(mtx, dispatch = TRUE)

			conditionals <- rJava::.jarray(Submatrix_temp[,1:(j-1)], dispatch = TRUE)
			
 			condEmbedDims2 <- rJava::.jarray(rep(condEmbedDims, j-1)) 
			cond_taus2 <- rJava::.jarray(rep(cond_taus, j-1))
			cond_delays2 <- rJava::.jarray(rep(cond_delays, j-1))
			conditionals2 <- rJava::.jarray(conditionals)

			collecTeCalc <- rJava::.jnew("mtinfodynamics/RunConditionalTransferEntropyCalculatorKraskov")
			rJava::.jcall(collecTeCalc, "V", "setProperty", "NOISE_LEVEL_TO_ADD", noise_level)
			rJava::.jcall(collecTeCalc,"V","initialise", k, k_tau, l, l_tau, delay, condEmbedDims2, cond_taus2, cond_delays2, num_permcheck)
			rJava::.jcall(collecTeCalc, "V", "runCTEKraskov", mtx_java, conditionals2, mode)

			for(t in 0:(num_timesteps-1)){
			  mtx <- result_mtx[[1]][[t+1]]
			  rmtx <- rJava::.jcall(collecTeCalc,"[[D",method = "getResultsMtxTimeStep", t)
			  vector_list_rmtx <- lapply(rmtx, function(mat) rJava::.jevalArray(mat, simplify = TRUE))
			  rmtx <- do.call(cbind, vector_list_rmtx)
			  mtx[i,j] <- rmtx[1,2]
			  result_mtx[[1]][[t+1]] <- mtx
			  mtx <- result_mtx[[2]][[t+1]]
			  Sgn <- rJava::.jcall(collecTeCalc,"[[D",method = "getSigMatx")
			  vector_list_sgn <- lapply(Sgn,function(mat) rJava::.jevalArray(mat, simplify = TRUE))
			  smtx <- do.call(rbind, vector_list_sgn)
			  mtx[i,j] <- smtx[1,2]
			  result_mtx[[2]][[t+1]] <- mtx
			}
			
      }

	 	}	

	}

	for(t in 1:num_timesteps){
	Adj <- result_mtx[[1]][[t]]
  rownames(Adj) <- rownames(data)
	result_mtx[[1]][[t]] <- Adj
	Sgn <- result_mtx[[2]][[t]]
	rownames(Sgn) <- rownames(data)
	result_mtx[[2]][[t]] <- Sgn
	  }

  	listnames <- colnames(data)
  	names(result_mtx[[1]]) <- listnames[seq(1,num_timesteps,1)]
  	names(result_mtx[[2]]) <- listnames[seq(1,num_timesteps,1)]
  	if(save){save(result_mtx,file = paste(filename,".Rdata", sep = ""))}
  	return(result_mtx)

	} else {

    DataSet <- .QtAC.split(num_timepoints,num_timesteps,Data1)      #Splitting & Extending

    # Calculation of adjacency and significance matrices:

    adjacency_matrices <- list()
    significance_matrices <- list()

    rJava::.jinit()
    rJava::.jaddClassPath(javapath)

    for (num_dataset in 1:length(DataSet)){

	      adj_matrix <- matrix(NaN, num_species, num_species-1)
	      sgn_matrix <- matrix(NaN, num_species, num_species-1)

	  for(i in 1:num_species){

    Submatrix <- DataSet[[num_dataset]]
		adjs <- rep(NaN, num_species-1)
		sgns <- rep(NaN, num_species-1)		

		destination <- Submatrix[,i]
		Submatrix_temp <- Submatrix[,-i]
		for(j in 1:(num_species-1)){

			source <- Submatrix_temp[,j]
			mtx <- cbind(source, destination)
			mtx_temp <- rJava::.jarray(mtx, dispatch = TRUE)

			if(j==1){
			  teCalc <- rJava::.jnew("mtinfodynamics/RunTransferEntropyCalculatorKraskov")
  			rJava::.jcall(teCalc,"V","setProperty", "NOISE_LEVEL_TO_ADD", noise_level)
  			rJava::.jcall(teCalc,"V","initialise",k,k_tau,l,l_tau,delay,num_permcheck)
  			rJava::.jcall(teCalc,"V","runTEKraskov",mtx_temp, "average")
  			Adj <- rJava::.jcall(teCalc,"[[D",method = "getResults")
  			Sgn <- rJava::.jcall(teCalc,"[[D",method = "getSigMatx")
  			vector_list_sgn <- lapply(Sgn,function(mat) rJava::.jevalArray(mat, simplify = TRUE))
  			vector_list_rmtx <- lapply(Adj,function(mat) rJava::.jevalArray(mat, simplify = TRUE))
  			smtx <- do.call(cbind, vector_list_sgn)
  			rmtx <- do.call(cbind, vector_list_rmtx)
  			adjs[j] <- rmtx[1,2]
			  sgns[j] <- smtx[1,2]
	   		}

			if(j>1){

			mtx <- cbind(source,destination)
			mtx_temp <- rJava::.jarray(mtx, dispatch = TRUE)

			conditionals2 <- rJava::.jarray(Submatrix_temp[,1:(j-1)], dispatch = TRUE)
			condEmbedDims2 <- rJava::.jarray(rep(condEmbedDims, j-1))
			cond_taus2 <- rJava::.jarray(rep(cond_taus, j-1))
			cond_delays2 <- rJava::.jarray(rep(cond_delays, j-1))

			collecTeCalc <- rJava::.jnew("mtinfodynamics/RunConditionalTransferEntropyCalculatorKraskov")
			rJava::.jcall(collecTeCalc, "V", "setProperty", "NOISE_LEVEL_TO_ADD", noise_level)
			rJava::.jcall(collecTeCalc,"V","initialise", k, k_tau, l, l_tau, delay, condEmbedDims2, cond_taus2, cond_delays2, num_permcheck)
			rJava::.jcall(collecTeCalc, "V", "runCTEKraskov", mtx_temp, conditionals2, "average")

			Adj <- rJava::.jcall(collecTeCalc,"[[D", method = "getResults")
			vector_list_rmtx <- lapply(Adj, function(mat) rJava::.jevalArray(mat, simplify = TRUE))
			rmtx <- do.call(cbind, vector_list_rmtx)

			Sgn <- rJava::.jcall(collecTeCalc,"[[D",method = "getSigMatx")
			vector_list_sgn <- lapply(Sgn,function(mat) rJava::.jevalArray(mat, simplify = TRUE))
			smtx <- do.call(cbind, vector_list_sgn)

			adjs[j] <- rmtx[1,2]
			sgns[j] <- smtx[1,2]
			
	  		}

	 	}

		adj_matrix[i,] <- adjs
		sgn_matrix[i,] <- sgns		

	}

  	rownames(adj_matrix) <- rownames(data)
	  rownames(sgn_matrix) <- rownames(data)

  	adjacency_matrices[[num_dataset]] <- adj_matrix
  	significance_matrices[[num_dataset]] <- sgn_matrix

  	Sys.time()     # just for time analysis

    }

    result_mtx[[1]] <- adjacency_matrices
    result_mtx[[2]] <- significance_matrices

    listnames <- colnames(data)
    names(result_mtx[[1]]) <- listnames[seq(num_timepoints,num_timesteps,1)]
    names(result_mtx[[2]]) <- listnames[seq(num_timepoints,num_timesteps,1)]
    if(save){save(result_mtx,file = paste(filename,".Rdata", sep = ""))}
    return(result_mtx)

	}	

}

#------------------------sign_CTE-----------------------------------
#' .signfactor_CTE
#'
#' This function sets all transfer entropies to 0 whose p-value is above the predefined significance level.
#' @export
#' @param result_mtx list of matrices containing summands of collective transfer entropies (one row per species) and corresponding significance matrices
#' @param signfac significance level
#' @return list of two lists (vectors of significant collective transfer entropies and significance matrices)
# @examples
# result_mtx_new <- QtAC.Signfactor(result_mtx,signfac = 0.05)

QtAC.signfactor_CTE <- function(result_mtx,signfac = 0.05){

  adjacent_matrices <- result_mtx[[1]]
  significance_matrices <- result_mtx[[2]]
  transformed_mtx <- list()
  adj_mtx <- list()

  for(num_mtx in 1:length(adjacent_matrices)){

    sign_matrix <- significance_matrices[[num_mtx]]
    adj_matrix <- adjacent_matrices[[num_mtx]]

    adj_matrix[sign_matrix > signfac] = 0

    adj_mtx[[num_mtx]] <- adj_matrix
  }
  
  vec <- list()
  
  for(i in 1:length(adjacent_matrices)){
    vec[[i]] <- rowSums(adj_mtx[[i]])
  }
  
  transformed_mtx[[1]] <- vec
  transformed_mtx[[2]] <- significance_matrices

  names(transformed_mtx[[1]]) <- names(result_mtx[[1]])
  names(transformed_mtx[[1]]) <- names(result_mtx[[2]])

  return(transformed_mtx)
}

  
#-------------------------AIS-------------------------------
#' QtAC_AIS 
#'
#' This function calculates the active information storage of each species for shifting time windows of fixed length or local values.
#' The output is a list of adjacency matrices and the corresponding significance matrices.
#' @export
#' @param data data array containing time series of the system's components' abundance data
#' @param num_timepoints length of the time windows of abundance data serving as basis of the active information storage estimations
#' @param javapath path of the file "MTinfodynamics.jar"
#' @param noise_level amount of random Gaussian noise added in the estimation
#' @param num_permcheck number of surrogate samples to bootstrap to generate the distribution in the significance test
#' @param k embedding length of destination past history to consider
#' @param k_tau embedding delay
#' @param mode Active information storage is either averaged over time windows ("average") or local values are estimated ("local"). In the latter case, the value of num_timepoints is ignored.
#' @param save If save=TRUE, the output is saved in a file called filename.
#' @return list of two vectors (AIS of all species and corresponding significance values)
# @examples
# result_mtx <- QtAC(Data,num_timepoints=6,JavaPath = "D:/Users/max.mustermann/Desktop/infoDynamics.jar",num_PermCheck=500)

QtAC_AIS <- function(data,num_timepoints,javapath,noise_level = "1e-20",num_permcheck=1000L,k=1L,k_tau=1L,mode = c("average","local"),save = FALSE, filename = "result_QtAC"){

  Data1 <- t(data)

  result_mtx <- list()
  num_timesteps <- dim(Data1)[1]
  num_species <- dim(Data1)[2]
  names_sp <- colnames(Data1)

  if(num_timepoints < 5 || num_timepoints > num_timesteps){
    print('num_timepoints has to be between 5 and the length of timepoints in the data')

  } else if(mode=="local"){ 

	DataSet <- .QtAC.split(num_timesteps,num_timesteps,Data1)      #Splitting & Extending

   	# Calculation of AIS and significance values:

	rJava::.jinit()
	rJava::.jaddClassPath(infodyn_path)

	Submatrix <- DataSet[[1]]

	mtx_java <- rJava::.jarray(Submatrix, dispatch = TRUE)

	aisCalc <- rJava::.jnew("mtinfodynamics/RunActiveInfoStorageCalculatorKraskov")
	rJava::.jcall(aisCalc,"V","setProperty", "NOISE_LEVEL_TO_ADD", noise_level)
	rJava::.jcall(aisCalc,"V","initialise",k,k_tau,num_permcheck)
	rJava::.jcall(aisCalc,"V","runAISKraskov", mtx_java, mode)

	rmtx <- rJava::.jcall(aisCalc,"[[D",method = "getMtx_time_steps")
	vector_list_rmtx <- lapply(rmtx, function(mat) rJava::.jevalArray(mat, simplify = TRUE))
	rmtx <- do.call(cbind, vector_list_rmtx)

	mtx_list <- list()
	for(i in 1:nrow(Submatrix)){
	  mtx <- rmtx[i,]
	  names(mtx) <- names_sp
	  mtx_list[[i]] <- mtx
	}

	Sgn <- rJava::.jcall(aisCalc,"[D",method = "getSigaisvector")
	names(Sgn) <- names_sp

	smtx_list <- list()
	for(i in 1:nrow(Submatrix)){
  	smtx_list[[i]] <- Sgn
	}

	result_mtx <- list(mtx_list, smtx_list)
	listnames <- colnames(data)
	names(result_mtx[[1]]) <- listnames[seq(1,num_timesteps,1)]
	names(result_mtx[[2]]) <- listnames[seq(1,num_timesteps,1)]
	if(save){save(result_mtx,file = paste(filename,".Rdata", sep = ""))}
	return(result_mtx)
	
} else{
	DataSet <- .QtAC.split(num_timepoints,num_timesteps,Data1)      #Splitting & Extending

        # Calculation of AIS and significance values:

        adjacency_matrices <- list()
        significance_matrices <- list()

    for (num_dataset in 1:length(DataSet)){

	Submatrix <- DataSet[[num_dataset]]
      
	rJava::.jinit()
 	rJava::.jaddClassPath(javapath)
  mtx_java <- rJava::.jarray(Submatrix, dispatch = TRUE)

	aisCalc <- rJava::.jnew("mtinfodynamics/RunActiveInfoStorageCalculatorKraskov")
	rJava::.jcall(aisCalc,"V","setProperty", "NOISE_LEVEL_TO_ADD", noise_level)
	rJava::.jcall(aisCalc,"V","initialise",k,k_tau,num_permcheck)
	rJava::.jcall(aisCalc,"V","runAISKraskov", mtx_java, "average")

	Adj <- rJava::.jcall(aisCalc,"[D",method = "getAisvector")
	names(Adj) <- names_sp

	Sgn <- rJava::.jcall(aisCalc,"[D",method = "getSigaisvector")
	names(Sgn) <- names_sp

	adjacency_matrices[[num_dataset]] <- Adj
  significance_matrices[[num_dataset]] <- Sgn

        Sys.time()     # just for time analysis
    }
    result_mtx[[1]] <- adjacency_matrices
    result_mtx[[2]] <- significance_matrices

  listnames <- colnames(data)
  names(result_mtx[[1]]) <- listnames[seq(num_timepoints,num_timesteps,1)]
  names(result_mtx[[2]]) <- listnames[seq(num_timepoints,num_timesteps,1)]
  if(save){save(result_mtx,file = paste(filename,".Rdata", sep = ""))}
  return(result_mtx)
 }

}


#-----------------------QtAC---------------------------------
#' QtAC (Main function)
#'
#' This function calculates the transfer entropy between two species each for shifting time windows of fixed length or local values.
#' The output is a list of adjacency matrices and the corresponding significance matrices.
#' @export
#' @param data data array containing time series of the system's components' abundance data
#' @param num_timepoints length of the time windows of abundance data serving as basis of the transfer entropy estimations
#' @param javapath path of the file "MTinfodynamics.jar"
#' @param noise_level amount of random Gaussian noise added in the estimation
#' @param num_permcheck number of surrogate samples to bootstrap to generate the distribution in the significance test
#' @param k embedding length of destination past history to consider
#' @param k_tau embedding delay for the destination variable
#' @param l embedding length of source past history to consider
#' @param l_tau embedding delay for the source variable
#' @param delay time lag between last element of source and destination next value
#' @param mode Active information storage is either averaged over time windows ("average") or local values are estimated ("local"). In the latter case, the value of num_timepoints is ignored.
#' @param save If save=TRUE, the output is saved in a file called filename.
#' @return list of two lists (adjacency and corresponding significance matrices)
# @examples
# result_mtx <- QtAC(Data,num_timepoints=6,JavaPath = "D:/Users/max.mustermann/Desktop/infoDynamics.jar",num_PermCheck=500)

QtAC <- function(data,num_timepoints,javapath,noise_level = "1e-20",num_permcheck=1000L,k=1L,k_tau=1L,l=1L,l_tau=1L,delay=1L, mode = c("average","local"), save = FALSE, filename = "result_QtAC"){

  Data1 <- t(data)

  result_mtx <- list()
  num_timesteps <- dim(Data1)[1]
  num_species <- dim(Data1)[2]

  if(num_timepoints < 5 || num_timepoints > num_timesteps){
    print('num_timepoints has to be between 5 and the length of timepoints in the data')

  } else if(mode == "local"){

    DataSet <- .QtAC.split(num_timesteps,num_timesteps,Data1)      #Splitting & Extending

    # Calculation of adjacency and significance matrices:

    adjacency_matrices <- list()
    significance_matrices <- list()

    result_mtx <- .QtAC.Kraskov(DataSet[[1]],num_species,noise_level,num_permcheck,k,k_tau,l,l_tau,delay,mode = "local",names_sp = rownames(data),javapath)

    Sys.time() # just for time analysis
    
    if(num_timesteps<15){
      vec <- 2
      for(i in 3:length(result_mtx[[1]])){
        if(i%%3!=1){
          vec <- c(vec,i)
        }
      }
        adj <- result_mtx[[1]]
        result_mtx[[1]] <- adj[-vec]
        sign <- result_mtx[[2]]
        result_mtx[[2]] <- sign[-vec]
    }
    
    listnames <- colnames(data)
    names(result_mtx[[1]]) <- listnames[seq(1,num_timesteps,1)]
    names(result_mtx[[2]]) <- listnames[seq(1,num_timesteps,1)]
    if(save){save(result_mtx,file = paste(filename,".Rdata", sep = ""))}
    return(result_mtx)

	} else {
    DataSet <- .QtAC.split(num_timepoints,num_timesteps,Data1)      #Splitting & Extending

    # Calculation of adjacency and significance matrices:

    adjacency_matrices <- list()
    significance_matrices <- list()

    for (num_dataset in 1:length(DataSet)){
      adj_sign <- .QtAC.Kraskov(DataSet[[num_dataset]],num_species,noise_level,num_permcheck,k,k_tau,l,l_tau,delay, mode = "average", names_sp = rownames(data),javapath)
      adjacency_matrices[[num_dataset]] <- adj_sign[[1]]
      significance_matrices[[num_dataset]] <- adj_sign[[2]]

      Sys.time()     # just for time analysis
    }
    result_mtx[[1]] <- adjacency_matrices
    result_mtx[[2]] <- significance_matrices

    listnames <- colnames(data)
    names(result_mtx[[1]]) <- listnames[seq(num_timepoints,num_timesteps,1)]
    names(result_mtx[[2]]) <- listnames[seq(num_timepoints,num_timesteps,1)]
    if(save){save(result_mtx,file = paste(filename,".Rdata", sep = ""))}
    return(result_mtx)
    
    }

}

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

  n <- 3*num_timepoints-2
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

.QtAC.Kraskov <- function(Submatrix,num_species,noise_level,num_permcheck,k,k_tau,l,l_tau,delay,mode=c("average","local"),names_sp,javapath){

  if(mode=="average"){

  rJava::.jinit()
  rJava::.jaddClassPath(javapath)
  mtx_java <- rJava::.jarray(Submatrix, dispatch = TRUE)

  teCalc <- rJava::.jnew("mtinfodynamics/RunTransferEntropyCalculatorKraskov")
  rJava::.jcall(teCalc,"V","setProperty", "NOISE_LEVEL_TO_ADD", noise_level)
  rJava::.jcall(teCalc,"V","initialise",k,k_tau,l,l_tau,delay,num_permcheck)
  rJava::.jcall(teCalc,"V","runTEKraskov", mtx_java, "average")

  rmtx <- rJava::.jcall(teCalc,"[[D",method = "getResults")
  vector_list_rmtx <- lapply(rmtx, function(mat) rJava::.jevalArray(mat, simplify = TRUE))
  Adj <- do.call(cbind, vector_list_rmtx)

  Adj[Adj<0] <- 0

  Sgn <- rJava::.jcall(teCalc,"[[D",method = "getSigMatx")
  vector_list_sgn <- lapply(Sgn,function(mat) rJava::.jevalArray(mat, simplify = TRUE))
  Sgn <- do.call(cbind, vector_list_sgn)

  colnames(Adj) <- names_sp
  rownames(Adj) <- names_sp
  colnames(Sgn) <- names_sp
  rownames(Sgn) <- names_sp

  adj_sign <- list(Adj,Sgn)

  return(adj_sign)
} else{

rJava::.jinit()
rJava::.jaddClassPath(javapath)

mtx_java <- rJava::.jarray(Submatrix, dispatch = TRUE)

teCalc <- rJava::.jnew("mtinfodynamics/RunTransferEntropyCalculatorKraskov")
rJava::.jcall(teCalc,"V","setProperty", "NOISE_LEVEL_TO_ADD", noise_level)
rJava::.jcall(teCalc,"V","initialise",k,k_tau,l,l_tau,delay,num_permcheck)
rJava::.jcall(teCalc,"V","runTEKraskov", mtx_java, "local")

t <- nrow(Submatrix) - 1
mtx_list <- list()
for(i in 0:t){
  rmtx <- rJava::.jcall(teCalc,"[[D",method = "getResultsMtxTimeStep", i)
  vector_list_rmtx <- lapply(rmtx, function(mat) rJava::.jevalArray(mat, simplify = TRUE))
  rmtx <- do.call(rbind, vector_list_rmtx)
  colnames(rmtx) <- names_sp
  rownames(rmtx) <- names_sp
  #rmtx[rmtx<0] <- 0
  mtx_list[[i+1]] <-rmtx
}

Sgn <- rJava::.jcall(teCalc,"[[D",method = "getSigMatx")
vector_list_sgn <- lapply(Sgn,function(mat) rJava::.jevalArray(mat, simplify = TRUE))
smtx <- do.call(cbind, vector_list_sgn)
colnames(smtx) <- names_sp
rownames(smtx) <- names_sp

smtx_list <- list()
for(i in 1:nrow(Submatrix)){
  smtx_list[[i]] <- smtx
}

  adj_sign <- list(mtx_list, smtx_list)

  return(adj_sign)

	}
}


#-----------------------sign------------------------------------
#' signfactor
#'
#' This function sets all entries in the adjacency matrices to 0 whose p-value is above the predefined significance level.
#' @export
#' @param result_mtx list of adjacency and significance matrices
#' @param signfac significance level
#' @return list of two lists (adjacency matrix containing only significant transfers and significance matrices)
# @examples
# result_mtx_new <- QtAC.Signfactor(result_mtx,signfac = 0.05)

QtAC.signfactor <- function(result_mtx,signfac = 0.05){

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

#---------------------------Potential---------------------------
# This function calculates the potential of a network.

.QtAC_potential <- function(mtx){

  mtx <- round(mtx,digits = 3)
  mtx[mtx<0] <- 0
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
  mtx[mtx<0] <- 0
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


.QtAC.resilience <- function(mtx, stand = "maxweight"){

  mtx <- round(mtx,digits = 3)
  mtx[mtx<0] <- 0

  length <- ncol(mtx)
  L_out <- pracma::zeros(length,length)
  L_in <- pracma::zeros(length,length)

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


  if(stand == "maxweight"){


    if(!(all(mtx==0))){

      b <- max(mtx)

      L_out <- (1/sqrt(b)) * L_out
      L_in <- (1/sqrt(b)) * L_in

    }

  }

  if(stand == "nodes"){


    if(!(all(mtx==0))){

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

  D_out <- round(abs(Re(eigen(L_out,only.values = TRUE)[[1]])),digits = 10)
  D_in <- round(abs(Re(eigen(L_in,only.values = TRUE)[[1]])),digits = 10)

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
#----------------------maturation--------------------------------------------
#' maturation
#'
#' This function computes the three systemic variables (potential, connectedness, and resilience) of each adjacency matrix.
#' @export
#' @param result_mtx list of adjacency matrices and significance matrices
#' @param res_stand standardization constant c of the Laplacian matrices ("none", "maxweight", "nodes", "maxweightnodes"). Let \eqn{N} be the number of nodes of the underlying graph and \eqn{M} its maximal edge weight. If res_stand = "none", \eqn{c = 1}. If res_stand = "maxweight", \eqn{c = \frac{1}{\sqrt{M}}}. If res_stand = "nodes", \eqn{c = \frac{\sqrt{N-1}}{N}}. If res_stand = "maxweightnodes", \eqn{c = \frac{\sqrt{N-1}}{N \cdot \sqrt{M}}}.
#' @return dataframe containing the three systemic variables (potential, connectedness, and resilience) of each adjacency matrix
# @examples
# Mat <- QtAC.maturation(result_mtx, res_stand = "maxweightnodes")

QtAC.maturation <- function(result_mtx, res_stand = "maxweight"){

  Adjacent_matrices <- result_mtx[[1]]
  length <- length(Adjacent_matrices)
  
  potential <- pracma::zeros(length,1)
  connectedness <- pracma::zeros(length,1)
  resilience <- pracma::zeros(length,1)

  for(num in 1:length){
    potential[num] <- .QtAC_potential(Adjacent_matrices[[num]])
    connectedness[num] <- .QtAC.connectedness(Adjacent_matrices[[num]])
    resilience[num] <- .QtAC.resilience(Adjacent_matrices[[num]], stand = res_stand)

  }

  maturation <- data.frame(potential,connectedness,resilience)
  rownames(maturation) <- names(result_mtx[[1]])
  return(maturation)

}

#----------------------network------------------------------------------
#' network
#'
#' This function plots a selected adjacency matrix as network.
#' @export
#' @param result_mtx list of adjacency matrices and significance matrices
#' @param num_mtx number of the adjacency matrix you want to plot (default = 1)
#' @param edge_label If edge_label = TRUE, the weight of the edges are plotted next to the edges.
#' @param dec number of decimal digits in the edge labels
#' @param layout layout format ("circle","star","fruchterman.reingold","grid","nicely")
#' @param edge_width muliplicator for the width of the edges
#' @param arrow_width muliplicator for the width of the arrows
#' @param col_node color of the vertices
#' @param col_edge color of the edges
#' @param vertex_label If vertex_label = "short", the first 4 letters of the components' names will be used as vertex labels, if vertex_label = "long", the whole names will be used as vertex labels. Via vertex_label =  c(...), customized names can be used as vertex labels.
#' @param title title of the network
#' @param save If save = TRUE, the network will be saved as a PNG file.
#' @param filename If save = TRUE, the network will be saved in a file called filename.
#' @return network plot and, if save = TRUE, a PNG file containing the plot
# @examples
# QtAC.network(result_mtx,3,edge_label = TRUE, dec = 3, layout = "nicely", edge_width = 30, arrow_width = 30,col_node = "red", vertex_label = "short", save = FALSE)

QtAC.network <- function(result_mtx,num_mtx = 1,edge_label = FALSE,dec = 2, layout = "nicely",edge_width = 3,arrow_width = 5,col_node = "palegreen3",col_edge = "steelblue3",vertex_label = "short", title = paste("Network of Adjacency matrix ",times[num_mtx]), save = FALSE,filename = paste("network_",num_mtx, sep = "")){

  adj_mtx <- result_mtx[[1]][[num_mtx]]
  
  adj_mtx[adj_mtx<0] <- 0

  graph_adj_mtx <- igraph::graph.adjacency(adj_mtx, mode="directed", weighted=TRUE)

  times <- names(result_mtx[[1]])

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

      vertex_label[name] <- paste(strsplit(names_long[name],"")[[1]][1],strsplit(names_long[name],"")[[1]][2],strsplit(names_long[name],"")[[1]][3],strsplit(names_long[name],"")[[1]][4],sep = "")

    }
  }

  #----------------------------------------------------

  if(save){
    filename <- paste(filename,'.png',sep='')
    png(filename,width = 800 ,height= 800,units = "px",type = "cairo")
    igraph::plot.igraph(graph_adj_mtx,layout=layout_graph, vertex.label = vertex_label, vertex.label.color = "black",vertex.color = col_node,edge.label = edge_label, edge.label.font= 2, edge.arrow.size =igraph::E(graph_adj_mtx)$weight/max(igraph::E(graph_adj_mtx)$weight)*arrow_width, edge.width=igraph::E(graph_adj_mtx)$weight/max(igraph::E(graph_adj_mtx)$weight)*edge_width, vertex.label.cex = 5, margin = 0, vertex.label.font=2, edge.curved= 0.1 , edge.color= col_edge, asp= 1)
    title(title, cex.main = 3)
    dev.off()
  }

  #plot the network, no matter whether it should be saved or not
  igraph::plot.igraph(graph_adj_mtx,layout=layout_graph, vertex.label = vertex_label, vertex.label.color = "black",vertex.color = col_node,edge.label = edge_label, edge.label.font= 2, edge.arrow.size =igraph::E(graph_adj_mtx)$weight/max(igraph::E(graph_adj_mtx)$weight)*arrow_width, edge.width=igraph::E(graph_adj_mtx)$weight/max(igraph::E(graph_adj_mtx)$weight)*edge_width, vertex.label.cex = 1, margin = 0, vertex.label.font=2, edge.curved= 0.1 , edge.color= col_edge, asp= 1)
  title(title, cex.main = 3)
}

#-------------------------2d---------------------------------------
#' 2dplot
#'
#' This function plots potential, connectedness and resilience with respect to time. Curves are interpolated via a piecewise cubic spline.
#' @export
#' @param mat dataframe containing a time series of systemic variables
#' @param prop If prop = NULL, the three systemic variables are plotted w.r.t time in one plot each. If prop = "potential", "connectedness", or "resilience", only the selected systemic variable is plotted w.r.t time.
#' @param time_int vector containing start time, end time, and step size to define the xaxis
#' @param time unit (i.e. "years", "steps",...)
#' @param save If save = TRUE, the 2D plot will be saved in a PNG file.
#' @param filename If save = TRUE, the network will be saved in a file called filename.
#' @return 2D plot and, if save = TRUE, a PNG file containing the plots.
# @examples
# QtAC.2dplot(maturation,prop = "potential", time_int = c(1990,2018,2))

QtAC.2dplot <- function(mat,prop = NULL, time_int = NULL, time = "time", save = FALSE, filename = "2dplot"){

  if(is.null(prop)){
    if(save){
      filename <- paste(filename,'.png',sep='')
      png(filename,width = 800 ,height= 800,units = "px",type = "cairo")
      #png(filename,units = "px",type = "cairo")
      par(mar=c(5,6,4,1)+.1)
      par(mfrow=c(3,1))
      .QtAC.potentialplot(mat[,1],time_int,rownames(mat),time,FALSE,"")
      .QtAC.connectednessplot(mat[,2],time_int,rownames(mat),time,FALSE,"")
      .QtAC.resilienceplot(mat[,3],time_int,rownames(mat),time,FALSE,"")
      dev.off()

      par(mar=c(5,6,4,1)+.1)
      par(mfrow=c(3,1))
      .QtAC.potentialplot(mat[,1],time_int,rownames(mat),time,FALSE,"")
      .QtAC.connectednessplot(mat[,2],time_int,rownames(mat),time,FALSE,"")
      .QtAC.resilienceplot(mat[,3],time_int,rownames(mat),time,FALSE,"")
    }
    else {
      par(mar=c(5,6,4,1)+.1)
      par(mfrow=c(3,1))
      .QtAC.potentialplot(mat[,1],time_int,rownames(mat), time,FALSE,"")
      .QtAC.connectednessplot(mat[,2],time_int,rownames(mat), time,FALSE,"")
      .QtAC.resilienceplot(mat[,3],time_int,rownames(mat), time,FALSE,"")
    }
  }

  else if(prop=="potential"){

    vec <- mat[,1]
    par(mar=c(5,6,4,1)+.1)
    par(mfrow = c(1,1))
    .QtAC.potentialplot(vec,time_int,rownames(mat),time,save,filename)

  }

  else if(prop=="connectedness"){

    vec <- mat[,2]
    par(mar=c(5,6,4,1)+.1)
    par(mfrow = c(1,1))
    .QtAC.connectednessplot(vec,time_int,rownames(mat),time,save,filename)
  }

  else if(prop=="resilience"){

    vec <- mat[,3]
    par(mar=c(5,6,4,1)+.1)
    par(mfrow = c(1,1))
    .QtAC.resilienceplot(vec,time_int,rownames(mat),time,save,filename)

  }

  else return("No available property")



}

#------------------------------Potentialplot------------------------------
# This function plots the potential w.r.t time. It is used in QtAC.2dplot. A curve is interpolated via a piecewise cubic spline.

.QtAC.potentialplot <- function(vec, time_int = NULL,names = seq(1,length(vec),1),time = "time",save,filename){


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

  if(save){
    filename <- paste(filename,'.png',sep='')
    png(filename,width = 800 ,height= 800,units = "px",type = "cairo")

    plot(xaxis_plot_spl,vec_ext,type = "l",col = "darkgreen",col.main = "darkblue",xlab = "", ylab = "", lwd = 3, axes = FALSE, frame.plot = TRUE)
    axis(side = 1, at = xaxis_values, labels = names, col.axis = "darkblue", cex.axis = 5)
    axis(side = 2, col.axis = "darkblue", cex.axis = 5)
    mtext(time,side = 1,line = 3, col = "darkblue", font = 2, cex = 2)
    mtext("potential",side = 2,line = 3, col = "darkgreen", font = 2, cex = 2)
    dev.off()
  }

  plot(xaxis_plot_spl,vec_ext,type = "l",col = "darkgreen",col.main = "darkblue",xlab = "", ylab = "", lwd = 3, axes = FALSE, frame.plot = TRUE)
  axis(side = 1, at = xaxis_values,labels = names, col.axis = "darkblue", cex.axis = 2)
  axis(side = 2, col.axis = "darkblue", cex.axis = 2)
  mtext(time,side = 1,line = 3, col = "darkblue", font = 2, cex = 2)
  mtext("potential",side = 2,line = 3, col = "darkgreen", font = 2, cex = 2)

}

#--------------------------Connectednessplot-------------------------
# This function plots the connectedness w.r.t time. It is used in QtAC.2dplot. A curve is interpolated via a piecewise cubic spline.

.QtAC.connectednessplot <- function(vec, time_int = NULL,names = seq(1,length(vec),1),time = "time",save,filename){


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


  if(save){
    filename <- paste(filename,'.png',sep='')
    png(filename,width = 800 ,height= 800,units = "px",type = "cairo")

    plot(xaxis_plot_spl,vec_ext,type = "l",col = "blue",col.main = "darkblue",xlab = "", ylab = "", lwd = 3, axes = FALSE, frame.plot = TRUE)
    axis(side = 1, at = xaxis_values,labels = names, col.axis = "darkblue", cex.axis = 2)
    axis(side = 2, col.axis = "darkblue", cex.axis = 2)
    mtext(time,side = 1,line = 3, col = "darkblue", font = 2, cex = 2)
    mtext("connectedness",side = 2,line = 3, col = "blue", font = 2, cex = 2)
    dev.off()
  }

  plot(xaxis_plot_spl,vec_ext,type = "l",col = "blue",col.main = "darkblue",xlab = "", ylab = "", lwd = 3, axes = FALSE, frame.plot = TRUE)
  axis(side = 1, at = xaxis_values,labels = names, col.axis = "darkblue", cex.axis = 2)
  axis(side = 2, col.axis = "darkblue", cex.axis = 2)
  mtext(time,side = 1,line = 3, col = "darkblue", font = 2, cex = 2)
  mtext("connectedness",side = 2,line = 3, col = "blue", font = 2, cex = 2)

}

#-----------------------------Resilienceplot-----------------------------
# This function plots the resilience w.r.t time. It is used in QtAC.2dplot. A curve is interpolated via a piecewise cubic spline.

.QtAC.resilienceplot <- function(vec, time_int = NULL,names = seq(1,length(vec),1), time = "time",save,filename){


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

  if(save){
    filename <- paste(filename,'.png',sep='')
    png(filename,width = 800 ,height= 800,units = "px",type = "cairo")

    plot(xaxis_plot_spl,vec_ext,type = "l",col = "red",col.main = "darkblue",xlab = "", ylab = "", lwd = 3, axes = FALSE, frame.plot = TRUE)
    axis(side = 1, at = xaxis_values,labels = names, col.axis = "darkblue", cex.axis = 2)
    axis(side = 2, col.axis = "darkblue", cex.axis = 2)
    mtext(time,side = 1,line = 3, col = "darkblue", font = 2, cex = 2)
    mtext("resilience",side = 2,line = 3, col = "red", font = 2, cex = 2)
    dev.off()
  }

  plot(xaxis_plot_spl,vec_ext,type = "l",col = "red",col.main = "darkblue",xlab = "", ylab = "", lwd = 3, axes = FALSE, frame.plot = TRUE)
  axis(side = 1, at = xaxis_values,labels = names, col.axis = "darkblue", cex.axis = 2)
  axis(side = 2, col.axis = "darkblue", cex.axis = 2)
  mtext(time,side = 1,line = 3, col = "darkblue", font = 2, cex = 2)
  mtext("resilience",side = 2,line = 3, col = "red", font = 2, cex = 2)

}

#----------------------2dmix---------------------------------------------
#' 2dmixplot
#'
#' This function plots two selected systemic variables w.r.t. each other. A curve is interpolated via a piecewise cubic spline.
#' @export
#' @param Mat data frame containing a time series of systemic variables
#' @param prop1 variable on x-axis ("potential","connectedness","resilience")
#' @param prop2 variable on y-axis ("potential","connectedness","resilience")
#' @param save If Save = TRUE, the 2D plot will be saved in a PNG file.
#' @param filename If Save = TRUE, the network will be saved in a file called filename.
#' @return 2D plot and, if Save = TRUE, a PNG file containing the plot
# @examples
# QtAC_2dmixplot(mat,"connectedness","resilience",save = TRUE, filename = "2dmixplot")

QtAC.2dmixplot <- function(mat,prop1, prop2,save = FALSE,filename = paste("2dmixplot_", prop1,"_", prop2)){

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
    return("Property 1 is not available")
  }

  vec1 <- mat[,prop_num1]
  vec2 <- mat[,prop_num2]

  vec1_ext <- pracma::pchip(0:(length(vec1)-1),vec1,seq(0,length(vec1)-1,by = 0.1))
  vec2_ext <- pracma::pchip(0:(length(vec2)-1),vec2,seq(0,length(vec2)-1,by = 0.1))

  filename <- paste(filename,".png",sep = "")

  par(mfrow = c(1,1))

  if(save){
    png(filename,width = 800 ,height= 800,units = "px",type = "cairo")

    plot(vec1_ext,vec2_ext,type = "l",col = color,main = paste(prop1," vs. ", prop2),xlab = prop1, ylab = prop2, lwd = 3)
    text(vec1,vec2,rownames(mat))
    dev.off()
  }

  plot(vec1_ext,vec2_ext,type = "l",col = color,main = paste(prop1," vs. ", prop2),xlab = prop1, ylab = prop2, lwd = 3)
  text(vec1,vec2,rownames(mat))
}


#--------------------3d-----------------------------------------
#' 3dplot
#'
#' 3D plot of the three systemic variables w.r.t each other. A curve is interpolated via a piecewise cubic spline.
#' @export
#' @param Mat data frame containing a time series of the three systemic variables
#' @param mat_points If mat_points = TRUE, the maturation points are visible.
#' @return 3D plot
# @examples
# QtAC.3dplot(mat = maturation, mat_points = TRUE)

QtAC.3dplot <- function(mat,mat_points=FALSE){

  axis <- seq(1,dim(mat)[1],1)
  steps <- seq(1,dim(mat)[1],0.01)

  x <- mat[,1]
  y <- mat[,2]
  z <- mat[,3]

  xspline <- pracma::pchip(axis,x,steps)
  yspline <- pracma::pchip(axis,y,steps)
  zspline <- pracma::pchip(axis,z,steps)

  #  x <- smooth.spline( x, spar=0 )
  #  y <- smooth.spline( y, spar=0 )
  #  z <- smooth.spline( z, spar=0 )

  #  x <- predict(x,steps)$y
  #  y <- predict(y,steps)$y
  #  z <- predict(z,steps)$y

  if(mat_points==TRUE){

    pointsize <- 0.5
  }
  else pointsize <- 0

  # rgl::plot3d(mat[,1], mat[,3], mat[,2], pch=19, cex=0.25, size=pointsize,type = "s", col="black",xlab = "potential",ylab = "resilience",zlab = "connectedness")
  rgl::plot3d(mat[,1], mat[,3], mat[,2], pch=19, cex=0.25, size=pointsize,type = "s", col="black",xlab = "",ylab = "",zlab = "")
  rgl::text3d(mat[,1],mat[,3], mat[,2], rownames(mat))

  rgl::lines3d(xspline,zspline,yspline,col = "orange",lwd=6)

}

