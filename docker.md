# Docker installation for running the QtAC container (optional)

**Docker installation in Debian based distributions (Ubuntu):**<br>
```console
$ sudo apt update

$ sudo apt install ca-certificates curl gnupg lsb-release

$ curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o \
     /usr/share/keyrings/docker-archive-keyring.gpg
$ echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] \
     https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

$ sudo apt-get update

$ sudo apt-get install docker-ce docker-ce-cli containerd.io
```

**Docker installation in Windows 10:**<br>
To install Docker in Windows 10 follow the steps [here](https://docs.docker.com/docker-for-windows/install/)<br>

**Docker installation in macOS:**<br>
To install Docker in macOS follow the steps [here](https://docs.docker.com/docker-for-mac/install/)<br>

**Test Docker installation:**<br>
Run this example in a terminal to make sure Docker is installed correctly:<br>
`docker run hello-world`

**NOTE**<br>
In linux you have to add sudo before docker!

----

# Run QtAC container and example

### Linux (Ubuntu):
* Download one time the image from Docker Hub. In a terminal copy and paste the following line:<br>
`docker pull hannahschrenk/qtac:latest`

* Run the container: The downloaded image contains a example dataset to test the features of the package.<br>
```console
sudo docker run -it -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix \
                                     --user="$(id --user):$(id --group)" \
                                     --device=/dev/dri:/dev/dri \
                                     -v $PWD:/output/ hannahschrenk/qtac:latest
```
The container will pop an R terminal where we can start using the QtAC package.

**NOTE**: <span style="color:red;">$PWD</span> is the current directory in the terminal (linux/macOS).<br>

### Windows 10/11:
* Download one time the image from Docker Hub. In a terminal copy and paste the following line:<br>
`docker pull hannahschrenk/qtac:latest`

### macOS (XX.XX):
* Download one time the image from Docker Hub. In a terminal copy and paste the following line:<br>
`docker pull hannahschrenk/qtac:latest`

Copy and paste.<br>

   1. INITIALIZATION:<br>
```R
library(QtAC)

work_folder  <- "/output/"                  # internal folder in the container mounted on $PWD
observ_data  <- "/Example/QtAC_SIHUMI.txt"  # internal file containing the example data in the container
infodyn_path <- "/dist/MTinfodynamics.jar"  # internal path of the MTinfodynamics.jar file in the container
setwd(work_folder)

num_timepoints <- 30   # length of the time windows serving as basis for the transfer entropy calculations
signfac <- 0.1         # significance level
```
   2. CALCULATIONS:<br>
```R
# load the data
Data <- QtAC.TXT.reader(observ_data, col_names = FALSE, row_names = TRUE)

# compute networks of information transfer for every time point starting from num_timepoints
result_mtx <- QtAC(Data, num_timepoints, infodyn_path, l = 10L, k = 10L, delay = 2L)

# take only information transfers passing the significance level into account
result_mtx_sig <- QtAC.signfactor(result_mtx, signfac)

# calculate the three systemic variables for every network
maturation <- QtAC.maturation(result_mtx_sig)
```
   3. VISUALIZATIONS:<br>
```R
# plot the first network of information transfers (corresponding to time point 30) and save it
QtAC.network(result_mtx_sig, num_mtx = 1, edge_label = TRUE, arrow_width = 2, layout = "nicely", save = TRUE)

# plot the development of potential, connectedness, and resilience over time and save it
QtAC.2dplot(maturation, save = TRUE)

# plot the development of potential and connectedness w.r.t. each other and save it
QtAC.2dmixplot(maturation, "potential", "connectedness", save = TRUE)

# plot a 3D plot of potential, connectedness, and resilience
QtAC.3dplot(maturation, mat_points = TRUE)
```
**NOTE**: All image files are saved in the current directory in the terminal.<br>

   4. EXIT CONTAINER:<br>
```R
quit()
```
