#' @export

simulationToDaily = function(times, states, dailyTimes, relevantStateNames)
{
	dailyStates = matrix(NA, length(dailyTimes), length(relevantStateNames))
	for(i in 1:length(dailyTimes))
	{
			ind = which(times <= dailyTimes[i])
			ind = length(ind)
			dailyStates[i, ] = states[ind, relevantStateNames]
	}
	return(dailyStates)
}