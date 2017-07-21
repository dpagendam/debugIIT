## debugIIT

## Stochastic modelling of Incompatible Insect Technique (IIT) programs for control of mosquito populations.
**Authors**: Dan Pagendam
**Contributors**: Nigel Snoad, Nigel Beebe, Brendan Trewin

debugIIT provides a set of functions for stochastic modelling of population dynamics using continuous-time Markov chains (birth-death processes).  The package is useful for stochastic simulation of possible population trajectories, which when aggregated, can be used to construct probabilistic predictive distributions for IIT program outcomes such as the likelihood of driving a population to extinction. The simulation code is largely written in C++ and makes use of the Rcpp package.

The package was developed at [CSIRO](http://www.csiro.au), Australia as part of collaboration with Verily's [Debug Project](https://debugproject.com/), targeting the mosquito _Aedes aegypti_.


### Package installation
First, clone the package from this repository using 

``` git clone https://github.com/dpagendam/debugIIT 
	cd debugIIT
	
```

You will be prompted to enter your GitHub username and password to complete the clone.

To install the package from GitHub, you will first need to install the devtools package in R using the command:

```install.packages("devtools")```

Once installed, you will need to load the devtools R package and install the friedsIndex R package using:

```
library(devtools)
devtools::install()
```

### Using this package

To use friedsIndex with some of the packaged example data, try:

```
library(adMRR)
library(raster)
# Load the example data for Innisfail, Queensland, Australia
data(landClass)
data(dataMRR6)
 
# Define the study region in landClass using the width in metres, the (lat, long) of 
# the upper left (UL) and lower right (LR) corners
width = 769.45
UL_lat = -1*(17 + 31/60 + 58.53/3600)
UL_long = 146 + 1/60 + 47.04/3600
LR_lat = -1*(17 + 32/60 + 29.53/3600)
LR_long = 146 + 2/60 + 13.18/3600
corners = matrix(rbind(c(UL_lat, UL_long), c(LR_lat, LR_long)),2, 2)
colnames(corners) = c("Lat", "Long")
corners_m = coordsToMetres(corners, UL_lat, UL_long, LR_lat, LR_long, width)
UL_m = corners_m[1, ]
LR_m = corners_m[2, ]
# Get the dimensions of the pixels in the study region
d = pixelDim(nrow(landClass), ncol(landClass), UL_m, LR_m)
dx = d["dx"]
dy = d["dy"]

# Add trap x, y locations to dataMRR1 data matrix
trapXY = coordsToMetres(dataMRR6[, c("Lat", "Long")], UL_lat, UL_long, LR_lat, LR_long, width)
dataMRR6 = cbind(dataMRR6, trapXY)

# Create a (lat, long) for the insect release site
releaseSite = matrix(c(-1*(17 + 32/60 + 18.77/3600), 146 + 1/60 + 58.16/3600), 1, 2, byrow = TRUE)

#releaseSite = matrix(NA, 5, 2)
#releaseSite[1, ] = c(-17.535362, 146.032838)
#releaseSite[2, ] = c(-17.536294, 146.032814)
#releaseSite[3, ] = c(-17.537180, 146.032790)
#releaseSite[4, ] = c(-17.538103, 146.032766)
#releaseSite[5, ] = c(-17.538806, 146.032748)

colnames(releaseSite) = c("Lat", "Long")

numInsectsReleased = 1250
#numInsectsReleased = rep(248, 5)
releaseSite_m = coordsToMetres(releaseSite, UL_lat, UL_long, LR_lat, LR_long, width)

# Create the spatial distribution of mosquito concentration (mosquitoes per square metre) for each pixel at the time of release.
initialConcentration = mrrInitialConditions(releaseSite, numInsectsReleased, nrow(landClass), ncol(landClass), UL_lat, UL_long, LR_lat, LR_long, dx, dy)

# Create a vector called times that contains the times at which the mosquito concentration grids (from the advection-diffusion model) are required.
deltaT = 1/24/4 #measured in days but equal to 15 minutes
tMax = 7 #measured in days
times = seq(0, tMax, deltaT)
 

# Construct Parameter Group for an advection-diffusion model
repellorLandClasses = c(1, 2, 4)
attractorLandClasses = c(5)
maxAdvection = 5.376
advectionDecayHalfDistance = 4.22
diffusivityValue = 68.415^2
deathRate = 0.382

diffusivity =  mrrParameter(gridType = "D", parameterNames = "diffusivity", lower = 0, upper = 40000, directions = "xy", relevantLandValues = c(1,2,3,4,5), parameterIsPixelValue = TRUE, paramValues = diffusivityValue)
attractors = mrrParameter(gridType = "v", parameterNames = c("maxAdvection", "advectionDecayHalfDistance"), lower = c(0, 0), upper = c(5, Inf), directions = "xy", relevantLandValues = attractorLandClasses, parameterIsPixelValue = FALSE, attractor = TRUE, paramValues = c(maxAdvection, advectionDecayHalfDistance))
repellors = mrrParameter(gridType = "v", parameterNames = c("maxAdvection", "advectionDecayHalfDistance"), lower = c(0, 0), upper = c(5, Inf), directions = "xy", relevantLandValues = repellorLandClasses, parameterIsPixelValue = FALSE, attractor = FALSE, paramValues = c(maxAdvection, advectionDecayHalfDistance))
K = mrrParameter(gridType = NULL, parameterNames = "K", lower = 0, upper = 1, directions = NA, relevantLandValues = NULL, parameterIsPixelValue = FALSE, attractor = FALSE, paramValues = 1)
bw = mrrParameter(gridType = NULL, parameterNames = "trapBandwidth", lower = 0 , upper = 25, directions = NA, relevantLandValues = NULL, parameterIsPixelValue = FALSE, attractor = FALSE, paramValues = NA)
mu = mrrParameter(gridType = NULL, parameterNames = "deathRate", lower = 0 , upper = 1.0, directions = NA, relevantLandValues = NULL, parameterIsPixelValue = FALSE, attractor = FALSE, paramValues = deathRate)
mrrParams = mrrParameterGroup(diffusivity, attractors, repellors, K, bw, mu)

# Run the model and store the output
run = adRun(mrrParams, initialConcentration, times, landClass, dx, dy, LR_m, UL_m, lrw = 16000000)

# Load the background image and plot the concentrations over time
bg_path = system.file("extdata", "Aerial.png", package="adMRR")
bg = raster(bg_path)

par(mfrow = c(4,2))
plot.concentration(run, 0, bg, title = "Day 0")
plot.concentration(run, 1, bg, title = "Day 1")
plot.concentration(run, 2, bg, title = "Day 2")
plot.concentration(run, 3, bg, title = "Day 3")
plot.concentration(run, 4, bg, title = "Day 4")
plot.concentration(run, 5, bg, title = "Day 5")
plot.concentration(run, 5, bg, title = "Day 6")
plot.concentration(run, 5, bg, title = "Day 7")

#Plot the advection field for the model
plot.advection(run, scale = 0.5, headlength = 0.025, thinning = 5)
```
