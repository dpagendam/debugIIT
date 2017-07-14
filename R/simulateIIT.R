#' @title Simulate an IIT program trajectory as a birth-death process.
#' 
#' @description \code{simulateIIT} simualtes a birth-death process model for a Wolbachia IIT release and returns a (stochastic) simulated trajectory for the system.
#' @param \code{params} A named parameter vector containing the model parameters.  This must include the birth rate "lambda"; the death rates of males and females "mu_m" and "mu_f"; the parameters "gamma_shape" and "gamma_rate" that govern the gamma distribution for time taken for immatures to develop from egg to adult; "K_eq" which is the number of adults we expect at equilibrium; "N_max" which is the ceiling on the number of immatures that can exist in the population at any point in time; and a Fried's index (denoted by the prefix "c_") for each type of insect in the population (e.g. "c_Wld" is the Fried's index for Wildtype males).  
#' @param \code{Wld_m} is the number of wildtype males that are expected at time zero.
#' @param \code{Wld_f} is the number of wildtype females that are expected at time zero.
#' @param \code{stochasticInitial} should the initial numbers be sampled treating \code{propTypes} as probabilities.
#' @param \code{numReleased} is the number of insects to be released at each each release time.  This should be a numeric vector of length equal to the number of release times.
#' @param \code{releaseMixture} a numeric vector containing the proportions of each insect type in the released mixture.  The vector should be named according to the insect types.
#' @param \code{contaminationProbs} a numeric vector containing the probabilities that the a released insect is a female rather than a male for each insect type released.
#' @param \code{propTypes} is the proportions of each type in the initial population.  This should be a numeric vector with names correcponding to the types of insects.
#' @param \code{releaseTimes} is a numeric vector containing the times at which the insects are released. 
#' @param \code{maxTime} is the maximum time to run the simulation for.  Simulations where the entire population goes extinct will be terminated earlier.
#' @param \code{maxSize} provides control over the blocksize used to extend the matrices used to store simulation outputs.
#' @return The function returns a named list containing a numeric vector of times ("t") and a matrix of corresponding states ("states").
#' @export


simulateIIT = function(params, Wld_m, Wld_f, stochasticInitial = TRUE, numReleased, releaseTypes, releaseMixture, contaminationProb, propTypes, releaseTimes, maxTime, maxSize = 1000000)
{
	params["theta"] = params["mu_m"]/params["mu_f"]
	times = c(releaseTimes, maxTime)
	n = round(params["gamma_shape"])
	types = names(propTypes)
	state = createBlankStartStateVector(params, types)
	names = names(propTypes)
	releasedTypes = names(releaseMixture)
	
	totalMales = 0;
	maleNumbers = c()
	maleNames = c()
	if(!stochasticInitial)
	{
		for(i in 1:length(names))
		{
			state[paste0(names[i], "_m")] = round(propTypes[i]*Wld_m)
			totalMales = totalMales + round(propTypes[i]*Wld_m)
		}
		
		for(i in 1:length(names))
		{
			for(j in 1:length(names))
			{
				state[paste0(names[i], "_f_", names[j])] = round(propTypes[i]*Wld_f*state[paste0(names[j], "_m")]/totalMales)
			}
		}
		
		for(i in 1:n)
		{
			for(j in 1:length(names))
			{
				state[paste0(names[j], "_imm_", i)] = round(propTypes[j]*params["K_eq"]*(params["mu_f"]*params["theta"] + params["mu_m"])/(params["gamma_rate"]*(1 + params["theta"])))
			}
		}
	}
	else
	{
		maleNumbers = rmultinom(1, round(Wld_m), propTypes)
		for(i in 1:length(names))
		{
			state[paste0(names[i], "_m")] = maleNumbers[i]
			totalMales = totalMales + maleNumbers[i]
		}
		
		femaleNumbers = rmultinom(1, round(Wld_f), propTypes)
		matingProbs = state[paste0(names, "_m")]/totalMales
		for(i in 1:length(names))
		{
			matingNumbers = rmultinom(1, femaleNumbers[i], matingProbs)
			for(j in 1:length(names))
			{
				state[paste0(names[i], "_f_", names[j])] = matingNumbers[j]
			}
		}
		
		for(i in 1:n)
		{
			immNumbers = rmultinom(1, round(params["K_eq"]*(params["mu_f"]*params["theta"] + params["mu_m"])/(params["gamma_rate"]*(1 + params["theta"]))), propTypes)
			for(j in 1:length(names))
			{
				state[paste0(names[j], "_imm_", i)] = immNumbers[j]
			}
		}
	}
	
	
	totalMalesReleased = 0
	numMalesByType = c()
	numReleasedByType = c()
	for(i in 1:length(releaseTypes))
	{
		numReleasedOfThisType = round(numReleased[1]*releaseMixture[releaseTypes[i]])
		malesReleasedOfThisType = rbinom(1, numReleasedOfThisType, 1 - contaminationProb)
		totalMalesReleased = totalMalesReleased + malesReleasedOfThisType
		numMalesByType = c(numMalesByType, malesReleasedOfThisType)
		numReleasedByType = c(numReleasedByType, numReleasedOfThisType)
		state[paste0(releaseTypes[i], "_m")] = state[paste0(releaseTypes[i], "_m")] + malesReleasedOfThisType
	}
	
	for(i in 1:length(releaseTypes))
	{
		femalesReleasedOfThisType = numReleasedByType[i] - numMalesByType[i]
		if(femalesReleasedOfThisType > 0)
		{
			femalesMatedByTypes = rmultinom(1, femalesReleasedOfThisType, maleNumbers/totalMales)
		}
		else
		{
			femalesMatedByTypes = rep(0, length(maleNumbers))
		}
		for(j in 1:length(releaseTypes))
		{
			state[paste0(releaseTypes[i], "_f_", releaseTypes[j])] = state[paste0(releaseTypes[i], "_f_", releaseTypes[j])] + femalesMatedByTypes[j]
		}
	}
	
	thisRound = simulateCTMC_cpp(state, types, params, times[1], times[2], maxSize, TRUE)
	
	t = thisRound$times
	state = thisRound$states

	immatureStateNames = c()
	for(i in 1:length(types))
	{
		immatureStateNames = c(immatureStateNames, paste0(types[i], "_imm_", 1:n))
	}
	femaleStateNames = c()
	for(i in 1:length(types))
	{
		for(j in 1:length(types))
		{
			femaleStateNames = c(femaleStateNames, paste0(types[i], "_f_", types[j]))
		}
	}
	colnames(state) <- c(paste0(types, "_m"), femaleStateNames, immatureStateNames)
	
	
	if(length(times) > 2)
	{
		for(x in 2:(length(times) - 1))
		{
			nextStart = state[nrow(state), ]
			
			totalMalesReleased = 0
			numMalesByType = c()
			numReleasedByType = c()
			for(i in 1:length(releaseTypes))
			{
				numReleasedOfThisType = round(numReleased[i]*releaseMixture[releaseTypes[i]])
				malesReleasedOfThisType = rbinom(1, numReleasedOfThisType, 1 - contaminationProb)
				totalMalesReleased = totalMalesReleased + malesReleasedOfThisType
				numMalesByType = c(numMalesByType, malesReleasedOfThisType)
				numReleasedByType = c(numReleasedByType, numReleasedOfThisType)
				nextStart[paste0(releaseTypes[i], "_m")] = nextStart[paste0(releaseTypes[i], "_m")] + malesReleasedOfThisType
			}
	
			for(i in 1:length(releaseTypes))
			{
				femalesReleasedOfThisType = numReleasedByType[i] - numMalesByType[i]
				if(femalesReleasedOfThisType > 0)
				{
					femalesMatedByTypes = rmultinom(1, femalesReleasedOfThisType, maleNumbers/totalMales)
				}
				else
				{
					femalesMatedByTypes = rep(0, length(maleNumbers))
				}
				
				for(j in 1:length(releaseTypes))
				{
					nextStart[paste0(releaseTypes[i], "_f_", releaseTypes[j])] = nextStart[paste0(releaseTypes[i], "_f_", releaseTypes[j])] + femalesMatedByTypes[j]
				}
			}	
			
			thisRound = simulateCTMC_cpp(nextStart, types, params, times[x], times[x + 1], maxSize, TRUE)
			t = c(t, thisRound$times)
			state = rbind(state, thisRound$states)
		}
	}
	
	return(list(t = t, state = state))
}

