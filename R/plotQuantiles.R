#' @title A function for visualising an object of type \code{trajectoryQuantileList}.
#' 
#' @description \code{fade} provides the name of a semi-transparent R colour after passing the name of a base colour or colours (e.g. "blue") and the degree of opacity desired.
#' @param \code{quantileList} is an object of class trajectoryQuantileList created using the function called trajectoryQuantiles. 
#' @param \code{times} is a numberic vector containing the times that the trajectories in \code{trajectoryQuantileList} correspond to.
#' @param \code{name} is the name of the state variable to be plotted (e.g. "Wld_m").
#' @param \code{col} A colour (preferably without added transparency) to use for plotting the trajectories.  Transparency will be added in the function.
#' @return The function returns a visualisation of temporal variability in a trajectory highlighting the 2.5th, 25th, 50th, 75th and 97.5th percetiles.
#' @export

plotQuantiles = function(quantileList, times, name, col)
{
	if(class(quantileList) != "trajectoryQuantileList")
	{
		stop("Input quantileList must be of class trajectoryQuantileList.  Use function trajectoryQuantiles to create this object.")
	}
	x = c(times, rev(times), times[1])
	plot(times, quantileList$percentile97.5[, name], col = "white", xlab = "time", ylab = name)
	y = c(quantileList$percentile2.5[, name], rev(quantileList$percentile97.5[, name]), quantileList$percentile2.5[, name][1])
	polygon(x = x, y = y, col = fade(col, 50), border = NA)
	y = c(quantileList$percentile25[, name], rev(quantileList$percentile75[, name]), quantileList$percentile25[, name][1])
	polygon(x = x, y = y, col = fade(col, 50), border = NA)
	lines(times, quantileList$percentile50[, name], col = col)
}




