#' @title Create a blank (named) state vector containing zero values.
#' 
#' @description \code{createBlankStartStateVector} uses the information about the types of insects and the parameters to construct a blank state vector.  This can then be populated with the correct numbers of each type to initiate a run of a model.
#' @param \code{params} A named parameter vector.  Must contain a parameter called "gamma_shape" for the time taken (assumed Gamma distributed with rate parameter equal to 1) for development from an egg to adult.
#' @param \code{types} is a character vector containing the types of insects.  One of these should be named "Wld" for the wildtype. 
#' @return The function returns a named numeric vector for the model.
#' @export



createBlankStartStateVector = function(params, types)
{
	n = round(params["gamma_shape"])
	numMaleStates = length(types)
	numFemaleStates = length(types)*(length(types) + 1)
	numImmatureStates = length(types)*n
	state = rep(0, numMaleStates + numFemaleStates + numImmatureStates )
	immatureStateNames = c()
	for(i in 1:length(types))
	{
		immatureStateNames = c(immatureStateNames, paste0(types[i], "_imm_", 1:n))
	}
	femaleStateNames = c()
	for(i in 1:length(types))
	{
		femaleStateNames = c(femaleStateNames, paste0(types[i], "_f_Unmated"))
		for(j in 1:length(types))
		{
			femaleStateNames = c(femaleStateNames, paste0(types[i], "_f_", types[j]))
		}
	}
	names(state) <- c(paste0(types, "_m"), femaleStateNames, immatureStateNames)
	return(state)
}

