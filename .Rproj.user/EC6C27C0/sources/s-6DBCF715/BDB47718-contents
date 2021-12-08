#---------------------------------------------------------------------------------------------#
#                            === EXAMPLE SIHUMI ===
#---------------------------------------------------------------------------------------------#
# Author: Hannah Schrenk and Nico Schreiber and Carlos Garcia-Perez
# Version: 0.1.1
#---------------------------------------------------------------------------------------------#
# INITIALIZATION

library(QtAC)

# Only for Windows: uncomment the following line and replace X with your Java version if a Java related error occurs
# Sys.setenv(JAVA_HOME="C:/Program Files/Java/jdk-XX.X.X/")

# set working directory (for Windows "C:/path/...")

work_folder  <- "path/to/folder"
setwd(work_folder)

# set path to QtAC_SIHUMI.txt

observ_data  <- "path/to/QtAC_SIHUMI.txt"

# set path to MTinfodynamics.jar (included in folder "dist")

infodyn_path <- "path/to/dist/MTinfodynamics.jar"

# set size of the time windows serving as basis for the transfer entropy calculations

num_timepoints <- 30

# set significance level

signfac        <- 0.1

#---------------------------------------------------------------------------------------------#
# CALCULATIONS

# load the data into the workspace

Data <- QtAC.TXT.reader(observ_data, col_names=FALSE, row_names = TRUE)

# compute networks of information transfer for every time point starting from num_timepoints

result_mtx <- QtAC(Data,num_timepoints, javapath = infodyn_path, l = 10L, k = 10L, delay = 2L, noise_level = "1e-20")

# take only information transfers passing the significance level into account

result_mtx_sig <- QtAC.signfactor(result_mtx,signfac)

# calculate the three systemic variables for every network

maturation <- QtAC.maturation(result_mtx_sig)

#----------------------------------------------------------------------------------------------#
# VISUALIZATIONS

# plot the first network of information transfers (corresponding to time point 30) and save it

QtAC.network(result_mtx_sig, num_mtx = 1, edge_label = TRUE, layout = "nicely", save = TRUE)

# plot the development of potential, connectedness, and resilience over time and save it

QtAC.2dplot(maturation, save = TRUE)

# plot the development of potential and connectedness w.r.t. each other

QtAC.2dmixplot(maturation, "potential", "connectedness", save = TRUE)

# plot a 3D plot of potential, connectedness, and resilience

QtAC.3dplot(maturation, mat_points = TRUE)
