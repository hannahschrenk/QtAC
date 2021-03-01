#---------------------------------------------------------------------------------------------#
#                            === EXAMPLE SIHUMI ===
#---------------------------------------------------------------------------------------------#
# Author: Nico Schreiber and Hannah Schrenk
# Version: 1.0
#---------------------------------------------------------------------------------------------#
# INITIALIZATION

work_folder  <- "D/Users/..."                           # folder containing the abundance data
observ_data  <- "QtAC_SIHUMI.txt"                       # file containing the abundance data
infodyn_path <- " D/Users/.../dist/MTinfodynamics.jar"  # path of the MTinfodynamics.jar file (included in "dist")
setwd(work_folder)

num_timepoints <- 30   # length of the time windows serving as basis for the transfer entropy calculations
signfac        <- 0.1  # significance level

#---------------------------------------------------------------------------------------------#
# CALCULATIONS

# load the data into the workspace

Data <- QtAC.TXT.reader(observ_data,col_names=FALSE,row_names = TRUE)

# compute networks of information transfer for every time point starting from num_timepoints

result_mtx <- QtAC(Data,num_timepoints, javapath = infodyn_path, l = 10L, k = 10L, delay = 2L)

# take only information transfers passing the significance level into account

result_mtx_sig <- QtAC.signfactor(result_mtx,signfac)

# calculate the three systemic variables for every network

maturation <- QtAC.maturation(result_mtx_sig)

#----------------------------------------------------------------------------------------------#
# VISUALIZATIONS

# plot the first network of information transfers (corresponding to time point 30) and save it

QtAC.network(result_mtx_sig,num_mtx = 1,edge_label = TRUE, arrow_width = 2, layout = "nicely", save = TRUE)

# plot the development of potential, connectedness, and resilience over time and save it

QtAC.2dplot(maturation, save = TRUE)

# plot the development of potential and connectedness w.r.t. each other

QtAC.2dmixplot(maturation, "potential", "connectedness")

# plot a 3D plot of potential, connectedness, and resilience

QtAC.3dplot(maturation,mat_points = TRUE)
