#' @title Simulate an IIT program trajectory as a birth-death process.
#' 
#' @description \code{createBlankStartStateVector} uses the information about the types of insects and the parameters to construct a blank state vector.  This can then be populated with the correct numbers of each type to initiate a run of a model.
#' @param \code{params} A named parameter vector.  Must contain a parameter called "gamma_shape" for the time taken (assumed Gamma distributed with rate parameter equal to 1) for development from an egg to adult.
#' @param \code{types} is a character vector.... 
#' @return The function returns a named numeric vector for the model.
#' @export

simulateCTMC = function(state, types, params, startTime, endTime, maxSize, store = TRUE)
{
	fpath_single = system.file("extdata", "simulateCTMC_cpp.R", package="debugIIT")
	source(fpath_single)
	result = debugIIT::simulateCTMC_cpp(state, types, params, startTime, endTime, maxSize, store)
	return(result)
}
