## debugIIT

## Stochastic modelling of Incompatible Insect Technique (IIT) programs for control of mosquito populations.
**Authors**: Dan Pagendam (CSIRO)

**Contributors**: Nigel Snoad (Verily), Nigel Beebe (University of Queensland), Brendan Trewin (CSIRO)

debugIIT provides a set of functions for stochastic modelling of population dynamics using continuous-time Markov chains (birth-death processes).  The package is useful for stochastic simulation of possible population trajectories, which when aggregated, can be used to construct probabilistic predictive distributions for IIT program outcomes such as the likelihood of driving a population to extinction. The simulation code is largely written in C++ and makes use of the Rcpp package.

The package was developed at [CSIRO](http://www.csiro.au), Australia as part of collaboration with Verily's [Debug Project](https://debugproject.com/), targeting the mosquito _Aedes aegypti_.


### Package installation
First, clone the package from this repository using 

``` git clone https://github.com/dpagendam/debugIIT ```

You will be prompted to enter your GitHub username and password to complete the clone.

To install the package then type

```	
	R CMD build debugIIT
	R CMD INSTALL debugIIT_1.0.tar.gz
```



### Using this package

The code below provides a test script to generate a single simulation in R

```
require(debugIIT)

# Scenario specific variables
numberOfBlocks = 1
contaminationProb = 1/250000
numReleasedPerHouse = 50
maxDays = 365
releaseDays = cumsum(c(0, rep(c(4,3), 25), 3))
densityDependentMating = FALSE


# Male and female death rates based on survival in Muir and Kay (1998)
mu_f_Wld = runif(1, 0.0943, 0.151)
mu_m_Wld = runif(1, 0.356, 0.562)

# Wolbachia death Rate modifier
wMelDeathRateModifier = runif(1, 1.0, 1.8)
wAlbDeathRateModifier = runif(1, 1.0, 1.4)
mu_f_WMel = mu_f_Wld*wMelDeathRateModifier
mu_m_WMel = mu_m_Wld
mu_f_WAlb = mu_f_Wld*wAlbDeathRateModifier
mu_m_WAlb = mu_m_Wld

# Parameters governing the delay while viable juveniles develop
gamma_shape = sample(10:18, size = 1)
gamma_rate = 1

# Fried's index for competition
c_Wld = 1
c_WAlb = runif(1, 0.7, 1.0)
c_WMel = runif(1, 0.7, 1.0)

# proportion that is wildtype versus WMel
propWildtype = runif(1, 0.1, 0.3)

# proportion of wildtype females that are in the mated state at equilibrium
propMated = runif(1, 0.2, 0.8)

###### Numbers of houses for reference #######
# Mourilyan has 256
# South Johnstone has 237
# Goondi Bend has 232
# Downtown Innisfail has 1561
##############################################

numHouses = sum(sample(15:30, size = numberOfBlocks))

# alpha is the multiplier applied to the total immatures at equilibrium and is used to 
# calculate maximum number of future adults (as immatures) that can be supported at a house
alpha = runif(numHouses, 1.0, 10.0)

# number of mosquitoes per house
numMosquitoes = sample(5:15, size = numHouses, replace = TRUE)

#Normal adult carrying capacity
K_eq = sum(numMosquitoes)

###########################################
###########################################
#			START - DONT MODIFY
###########################################

# These are other parameters that are derived from equilibria
theta = mu_m_Wld/mu_f_Wld
Wld_m = K_eq/(1 + theta)
Wld_f_Unmated = (1 - propMated)*K_eq*theta/(1 + theta)
Wld_f_Wld = propMated*K_eq*theta/(1 + theta)
I_eq = K_eq*(mu_f_Wld*theta + mu_m_Wld)/(gamma_rate*(1 + theta))
I_max = sum(round(I_eq*gamma_shape/numHouses)*alpha)
lambda = I_max*gamma_rate*I_eq*(1 + theta)/(K_eq*theta*propMated*(I_max - I_eq*gamma_shape))
eta = numHouses/Wld_m*((gamma_rate*I_eq*(1 + theta))/(2*(1 - propMated)*K_eq*theta) - mu_f_Wld)

###########################################
###########################################
#			END - DONT MODIFY
###########################################


# total released for the mixed area for the IIT
numReleased = numReleasedPerHouse*numHouses

DD_mating = ifelse(densityDependentMating, 1, 0)

params = c(mu_f_Wld, mu_f_WMel, mu_f_WAlb, mu_m_Wld, mu_m_WMel, mu_m_WAlb, gamma_shape, gamma_rate, c_Wld, c_WMel, c_WAlb, I_max, lambda, K_eq, eta, numHouses, DD_mating)
names(params) <- c("mu_f_Wld", "mu_f_WMel", "mu_f_WAlb", "mu_m_Wld", "mu_m_WMel", "mu_m_WAlb", "gamma_shape", "gamma_rate", "c_Wld", "c_WMel", "c_WAlb", "I_max", "lambda", "K_eq", "eta" , "H", "DD_mating")

propTypes = c(propWildtype, (1 - propWildtype), 0)
names(propTypes) = c("Wld", "WMel", "WAlb")
releaseMixture = c(0, 0, 1)
names(releaseMixture) = c("Wld", "WMel", "WAlb")


###########################################
###########################################
#		SIMULATE SOME TRAJECTORIES
###########################################
```

In the example below, we're going to run the model as a metapopulation model with two blocks.

```
# The times and states for the simulations are stored in the lists below
tList = list()
stateList = list()

#Number of simulations
numSims = 1

suppressionObjective = 0.1*Wld_m
pastSuppressionModifier = 0.25
releaseNumbers = rep(numReleased, length(releaseDays))

# This matrix defines the rates at which adults transition between the two blocks.
# The rows and columns correspond to the two blocks and the rates are for movement from
# the row block to the column block.
blockImmigrationRates = matrix(0, 2, 2)
blockImmigrationRates[1, 2] = 0.25
blockImmigrationRates[2, 1] = 0.25

# You can use different parameters on each block, but we'll use the same ones for the two blocks below.
paramsList = list()
paramsList[[1]] = params
paramsList[[2]] = params

# We'll also specify the same initial numbers of wildtype male and females below
Wld_m = rep(Wld_m, 2)
Wld_f_Unmated = rep(Wld_f_Unmated, 2)
Wld_f_Wld= rep(Wld_f_Wld, 2)

# We'll also use the same release days for both blocks
numReleased = rep(numReleased, length(releaseDays))
numReleased = rbind(numReleased, numReleased)

simulation = simulateIIT_metapopulation(paramsList, blockImmigrationRates, Wld_m, Wld_f_Unmated, Wld_f_Wld, stochasticInitial = TRUE, numReleased = numReleased, ratioReleased = NULL, releaseMixture, contaminationProb, propTypes, releaseTimes = releaseDays, maxTime  = maxDays, initState = NULL, maxSize = 1000000)
# Create a plot for just the wildtype males for example
plot(simulation$t, simulation$state[, "Wld_m.Block_1"], ty= "l", xlab = "Days", ylab = "Number of Wildtype Males on Block 1")
```
If you'd now like to generate and visualise an ensemble of trajectories, you can try the code below:

```
# Create an ensemble of simualations that have been reduced to a daily timestep
numSims = 10
regularTimes = 0:365
relevantStateNames = c("Wld_m.Block_1", "WMel_m.Block_1", "WAlb_m.Block_1", "Wld_Unmated.Block_1", "Wld_f_Wld.Block_1", "Wld_f_WMel.Block_1", "Wld_f_WAlb.Block_1", "WMel_Unmated.Block_1", "WMel_f_Wld.Block_1", "WMel_f_WMel.Block_1", "WMel_f_WAlb.Block_1", "WAlb_Unmated.Block_1", "WAlb_f_Wld.Block_1", "WAlb_f_WMel.Block_1", "WAlb_f_WAlb.Block_1")
dailySimStorage = array(NA, c(numSims, length(regularTimes), length(relevantStateNames)))


for(i in 1 :numSims)
{
	cat("Performing simulation ", i, "....")
	simulation = simulateIIT(params = params, Wld_m = Wld_m, Wld_f_Unmated = Wld_f_Unmated, Wld_f_Wld = Wld_f_Wld, stochasticInitial = TRUE, numReleased = rep(numReleased, length(releaseDays)), ratioReleased = NULL, releaseMixture = releaseMixture, contaminationProb = contaminationProb, propTypes = propTypes, releaseTimes = releaseDays, maxTime = maxDays, maxSize = 1000000)
	dailySimStorage[i, , ] = simulationToDaily(simulation$t, simulation$state, regularTimes, relevantStateNames)
	cat("Done. \n")
}


quantileList = trajectoryQuantiles(dailySimStorage, relevantStateNames)

pdf("quantileTrajectories.pdf")
par(mfrow = c(3,1))
plotQuantiles(quantileList, regularTimes, "Wld_m.Block_1", "red")
plotQuantiles(quantileList, regularTimes, "WMel_m.Block_1", "purple")
plotQuantiles(quantileList, regularTimes, "WAlb_m.Block_1", "blue")

par(mfrow = c(3,1))
plotQuantiles(quantileList, regularTimes, "Wld_f_Wld.Block_1", "red")
plotQuantiles(quantileList, regularTimes, "Wld_f_WMel.Block_1", "purple")
plotQuantiles(quantileList, regularTimes, "Wld_f_WAlb.Block_1", "blue")

par(mfrow = c(3,1))
plotQuantiles(quantileList, regularTimes, "WMel_f_Wld.Block_1", "red")
plotQuantiles(quantileList, regularTimes, "WMel_f_WMel.Block_1", "purple")
plotQuantiles(quantileList, regularTimes, "WMel_f_WAlb.Block_1", "blue")

par(mfrow = c(3,1))
plotQuantiles(quantileList, regularTimes, "WAlb_f_Wld.Block_1", "red")
plotQuantiles(quantileList, regularTimes, "WAlb_f_WMel.Block_1", "purple")
plotQuantiles(quantileList, regularTimes, "WAlb_f_WAlb.Block_1", "blue")

dev.off()

```

