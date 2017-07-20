#' @export

trajectoryQuantiles = function(dailySim, relevantStateNames)
{
	d = dim(dailySim)
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
			percentile97.5[i, j] = quantile(as.numeric(dailySim[ , i, j]), 0.975)
			percentile75[i, j] = quantile(as.numeric(dailySim[ , i, j]), 0.75)
			percentile50[i, j] = quantile(as.numeric(dailySim[ , i, j]), 0.5)
			percentile25[i, j] = quantile(as.numeric(dailySim[ , i, j]), 0.25)
			percentile2.5[i, j] = quantile(as.numeric(dailySim[ , i, j]), 0.025)
		}
	}
	return(list(percentile97.5 = percentile97.5, percentile75 = percentile75, percentile50 = percentile50, percentile25 = percentile25, percentile2.5 = percentile2.5))
}


