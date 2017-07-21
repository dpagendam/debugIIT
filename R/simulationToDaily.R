#' @title Convert a simulated birth-death process trajectory into a multivariate timeseries at specified time points.
#' 
#' @description \code{simulationToDaily} creates a matrix whose rows are the state of a simulated birth-death process at specific points in time.  The function is useful for turning a simulation into a multivariate timeseries on a regular timestep for example.
#' @param \code{times} a numeric vector containing the event occurence times in a simulated birth-death process.
#' @param \code{states} a numeric matrix containing (as rows) the state vectors of the simulated birth-death process at the times given in \code{times}.  The column names of the matrix should contain the state names of each element in the state vector.
#' @param \code{newTimes} a numeric vector containing the times at which the constructed multivariate time series should be created (can be at regular or irregular intervals).
#' @param \code{relevantStateNames} a character vector containing the names of the state variables that should be included in the returned in the multivariate timeseries. These names should be elements in the set of column names in \code{states}.
#' @return a matrix object with rows corresponding to each of the time points in the input \code{times}.
#' @export

simulationToDaily = function(times, states, newTimes, relevantStateNames)
{
	dailyStates = matrix(NA, length(newTimes), length(relevantStateNames))
	for(i in 1:length(newTimes))
	{
			ind = which(times <= newTimes[i])
			ind = length(ind)
			dailyStates[i, ] = states[ind, relevantStateNames]
	}
	return(dailyStates)
}