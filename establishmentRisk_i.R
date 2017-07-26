require(debugIIT)
# sex separation of 1:250000
# single block
setwd("~/establishmentRisk/scenario10/")


args = commandArgs(trailingOnly = TRUE)
argVal = as.integer(args[1])
set.seed(argVal)

# Male and female death rates
mu_f = runif(1, 1/15, 1/3)
maleModifier = runif(1, 1.0, 2.0)
mu_m = maleModifier*mu_f

# Parameters governing the delay while viable juveniles develop
gamma_shape = sample(10:18, size = 1)
gamma_rate = 1

# Fried's index for competition
c_Wld = 1
c_WAlb = runif(1, 0.7, 1.0)
c_WMel = runif(1, 0.7, 1.0)

# proportion that is wildtype versus WMel
propWildtype = runif(1, 0.1, 0.3)

# Number of houses
# Mourilyan has 256
# South Johnstone has 237
# Goondi Bend has 232
# Downtown Innisfail has 1561
numHouses = sample(15:30, size = 10)

# number of mosquitoes per house
numMosquitoes = sample(5:15, size = sum(numHouses), replace = TRUE)

#Normal adult carrying capacity
K_eq = sum(numMosquitoes)

# 90% suppression level
suppressionObjective = K_eq*0.1
pastSuppressionModified = 0.25

###########################################
###########################################
#			OTHER - DONT MODIFY
###########################################

# Other
theta = mu_m/mu_f
Wld_m = K_eq/(1 + theta)
Wld_f = K_eq*theta/(1 + theta)
I_eq = K_eq*(mu_f*theta + mu_m)/(gamma_rate*(1 + theta))
N_multiplier = runif(1, 1.0, 10.0)
N_max = gamma_shape*I_eq*N_multiplier
lambda = N_max*gamma_rate*I_eq*(1 + theta)/(K_eq*theta*(N_max - I_eq*gamma_shape))

# overburden to use for the SIT
releaseDays = cumsum(c(0, rep(c(4,3), 25), 3))
maxDays = 365
numReleased = 50*numHouses

#contamination rate
contaminationProb = 1/250000

params = c(mu_f, mu_m, gamma_shape, gamma_rate, c_Wld, c_WMel, c_WAlb, N_max, lambda, K_eq)
names(params) <- c("mu_f", "mu_m", "gamma_shape", "gamma_rate", "c_Wld", "c_WMel", "c_WAlb", "N_max", "lambda", "K_eq")


###########################################
###########################################
#		SIMULATE SOME TRAJECTORIES
###########################################

#The times and states for the simulations are stored in the lists below
tList = list()
stateList = list()

#Number of simulations
numSims = 1


#Do the simulations using the Rcpp modules
propTypes = c(propWildtype, (1 - propWildtype), 0)
names(propTypes) = c("Wld", "WMel", "WAlb")
releaseMixture = c(0, 0, 1)
names(releaseMixture) = c("Wld", "WMel", "WAlb")
releaseNumbers = rep(numReleased, length(releaseDays))

simulation = list()
simulation1 = simulateIIT(params, Wld_m, Wld_f, TRUE, numReleased = rep(numReleased, length(releaseDays)), ratioReleased = NULL, releaseMixture, contaminationProb, propTypes, releaseDays, maxDays, initState = NULL, maxSize = 1000000)

malePop = simulation1$state[, "Wld_m"] + simulation1$state[, "WMel_m"]
objAchievedInd = which(malePop < suppressionObjective)
if(length(objAchievedInd) == 0)
{
	simulation = simulation1
}
if(length(objAchievedInd) > 0)
{
	objAchievedInd = objAchievedInd[1]
	objAchievedTime = simulation1$t[objAchievedInd]
	objectiveState = simulation1$state[objAchievedInd, ]
	remainingReleaseDays = c(0, releaseDays[releaseDays >= objAchievedTime] - objAchievedTime)
	remainingReleaseNumbers = round(pastSuppressionModified*releaseNumbers[releaseDays >= objAchievedTime])
	simulation2 = simulateIIT(params, Wld_m, Wld_f, TRUE, numReleased = remainingReleaseNumbers, ratioReleased = NULL, releaseMixture, contaminationProb, propTypes, remainingReleaseDays, maxDays - objAchievedTime, initState = objectiveState, maxSize = 1000000)
	simulation$t = c(simulation1$t[1:objAchievedInd], simulation2$t[-1] + objAchievedTime)
	simulation$state = rbind(simulation1$state[1:objAchievedInd, ], simulation2$state[-1, ])
}

save(list = c("simulation", "params", "contaminationProb", "propWildtype", "releaseDays", "Wld_m", "Wld_f", "numReleased", "numHouses", "K_eq", "numMosquitoes"), file = paste0("Simulation_", argVal, ".RData"))


#pdf("TrajectoryPlots.pdf")
#for(i in 1:45)
#{
#	plot(simulationList[[1]]$t, simulationList[[1]]$state[, i], ty = "l", ylab = colnames(simulationList[[1]]$state)[i], col = fade("black", 50), ylim = c(0, max(simulationList[[1]]$state[, i])*1.2))
#	for(j in 2:length(simulationList))
#	{
#		lines(simulationList[[j]]$t, simulationList[[j]]$state[, i], col = fade("black", 50))
#	}
#}
#dev.off()

