#' @export

plotQuantiles = function(quantileList, times, name, col)
{
	x = c(times, rev(times), 0)
	plot(times, quantileList$percentile97.5[, name], col = "white", xlab = "time", ylab = name)
	y = c(quantileList$percentile2.5[, name], rev(quantileList$percentile97.5[, name]), quantileList$percentile2.5[, name][0])
	polygon(x = x, y = y, col = fade(col, 50))
	y = c(quantileList$percentile25[, name], rev(quantileList$percentile75[, name]), quantileList$percentile25[, name][0])
	polygon(x = x, y = y, col = fade(col, 50))
	lines(0:365, quantileList$percentile50[, name], col = col)
}




