#' @title Create an object of class \code{trajectoryQuantileList} that summarises the variability in an ensemble of simulated trajectories.
#' 
#' @description \code{trajectoryQuantiles} returns an object of class \code{trajectoryQuantileList} that contains named elements that summarise, at each point in time, the 2.5th, 25th, 50th, 75th and 97.5th percentiles of an ensemble of trajectories.
#' @param \code{regSimArray} an array object that stores an ensemble of simulated trajectories that share a regular timestep.  The array must have three dimensions: the first corresponds to the number of simulations in the ensemble; the second corresponds to the number of time-points in the simulated trajectories (these must have been regularised to a common set of times); and the third being the number of states in the state vector of the simulated process.
#' @param \code{relevantStateNames} a character vector containing the names of the elements in the state vector for \code{regSimArray} (i.e. the state names).
#' @return a matrix object with rows corresponding to each of the time points in the input \code{times}.
#' @export

trajectoryQuantiles = function(regSimArray, relevantStateNames)
{
	d = dim(regSimArray)
	numSims = d[1]
	numDays = d[2]
	numState = d[3]
	
	percentile97.5 = matrix(NA, numDays, numState)
	percentile75 = matrix(NA, numDays, numState)
	percentile50 = matrix(NA, numDays, numState)
	percentile25 = matrix(NA, numDays, numState)
	percentile2.5 = matrix(NA, numDays, numState)
	
	colnames(percentile97.5) = relevantStateNames
	colnames(percentile75) = relevantStateNames
	colnames(percentile50) = relevantStateNames
	colnames(percentile25) = relevantStateNames
	colnames(percentile2.5) = relevantStateNames

	for(i in 1:numDays)
	{
		for(j in 1:numState)
		{
			percentile97.5[i, j] = quantile(as.numeric(regSimArray[ , i, j]), 0.975)
			percentile75[i, j] = quantile(as.numeric(regSimArray[ , i, j]), 0.75)
			percentile50[i, j] = quantile(as.numeric(regSimArray[ , i, j]), 0.5)
			percentile25[i, j] = quantile(as.numeric(regSimArray[ , i, j]), 0.25)
			percentile2.5[i, j] = quantile(as.numeric(regSimArray[ , i, j]), 0.025)
		}
	}
	result = list(percentile97.5 = percentile97.5, percentile75 = percentile75, percentile50 = percentile50, percentile25 = percentile25, percentile2.5 = percentile2.5)
	class(resut) = "trajectoryQuantileList"
	return()
}


