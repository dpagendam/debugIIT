require(debugIIT)


#args = commandArgs(trailingOnly = TRUE)
#argVal = as.integer(args[1])
#set.seed(argVal)

argVal = 1

# Scenario specific variables
numberOfBlocks = 1
contaminationProb = 1/10000
overflooding = 5
maxDays = 365*2
releaseDays = cumsum(c(0, rep(c(2, 2, 3), 52)))
densityDependentMating = FALSE;




goodParams = FALSE
while(!goodParams)
{
	#We use pripor distributions over the parameters, but not all combinations of these will be valid.
	# Some rejection sampling is used to ensure that the sampled parameters are valid.

	# Male and female death rates based on survival in Muir and Kay (1998)
	mu_f_Wld = runif(1, 0.0943, 0.151)
	mu_m_Wld = runif(1, 0.2231, 0.562)

	# Wolbachia death Rate modifier
	wMelDeathRateModifier = runif(1, 1.0, 1.0) #1.8 from Axford lab study
	wAlbDeathRateModifier = runif(1, 1.0, 1.0) #1.4 from Axford lab study
	mu_f_WMel = mu_f_Wld*wMelDeathRateModifier
	mu_m_WMel = mu_m_Wld
	mu_f_WAlb = mu_f_Wld*wAlbDeathRateModifier
	mu_m_WAlb = mu_m_Wld

	# Parameters governing the delay while viable juveniles develop
	gamma_shape = sample(10:40, size = 1)
	gamma_rate = 1

	# Fried's index for competition
	c_Wld = 1
	c_WAlb = runif(1, 0.7, 1.0)
	c_WMel = runif(1, 0.7, 1.0)

	# proportion that is wildtype
	propWildtype = 1

	# proportion of wildtype females that are in the mated state at equilibrium
	propMated = runif(1, 0.2, 0.8)

	#The intrinsic reproductive rate of production of males and females
	lambda = runif(1, 0.2, 0.6)
	
	
	###### Numbers of houses for reference #######
	# Mourilyan has 256
	# South Johnstone has 237
	# Goondi Bend has 232
	# Downtown Innisfail has 1561
	##############################################

	numHouses = 20

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

	I_max = gamma_shape*I_eq/(1 - (gamma_rate*I_eq)/(lambda*Wld_f_Wld))

	if(!densityDependentMating)
	{
		eta = ((gamma_rate*I_eq*(1 + theta))/(2*(1 - propMated)*K_eq*theta) - mu_f_Wld)
	}
	if(densityDependentMating)
	{
		eta = numHouses/Wld_m*((gamma_rate*I_eq*(1 + theta))/(2*(1 - propMated)*K_eq*theta) - mu_f_Wld)
	}
	###########################################
	###########################################
	#			END - DONT MODIFY
	###########################################
	
	Wld_m = round(Wld_m)
	Wld_f_Unmated = round(Wld_f_Unmated)
	Wld_f_Wld = round(Wld_f_Wld)
	I_eq = round(I_eq)
	I_max = round(I_max)
	
	#Check Parametric Constriants
	test1 = gamma_rate*I_eq/(lambda*Wld_f_Wld)
	test2_1 = 0.5*gamma_rate*I_eq
	test2_2 = mu_f_Wld*Wld_f_Unmated
	
	if(test1 < 1 & test2_1 > test2_2)
	{
		goodParams = TRUE
	}

}


# total released for the mixed area for the IIT
numReleased = Wld_m*overflooding

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


#The times and states for the simulations are stored in the lists below
tList = list()
stateList = list()

#Number of simulations
numSims = 1

thisSim = simulateIIT(params = params, Wld_m = Wld_m, Wld_f_Unmated = Wld_f_Unmated, Wld_f_Wld = Wld_f_Wld, stochasticInitial = TRUE, numReleased = NULL, ratioReleased = overflooding, releaseMixture = releaseMixture, contaminationProb = contaminationProb, propTypes = propTypes, releaseTimes = releaseDays, maxTime = maxDays, maxSize = 1000000)
simulation = list()
simulation$t = thisSim$t
simulation$state = thisSim$state


save(list = c("simulation", "params", "contaminationProb", "propWildtype", "releaseDays", "Wld_m", "Wld_f_Wld", "Wld_f_Unmated", "numReleased", "numHouses", "K_eq", "numMosquitoes"), file = paste0("Simulation_", argVal, ".RData"))


#pdf("TrajectoryPlots.pdf")
#for(i in 1:48)
#{
#	plot(simulation$t, simulation$state[, i], ty = "l", ylab = colnames(simulation$state)[i], col = fade("black", 50), ylim = c(0, max(simulation$state[, i])*1.2))
#}
#dev.off()

