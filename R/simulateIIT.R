#' @title Simulate an IIT program trajectory as a birth-death process.
#' 
#' @description \code{simulateIIT} simualtes....
#' @param \code{params} A named parameter vector.  
#' @return The function returns a ....
#' @export


simulateIIT = function(params, Wld_m, Wld_f, numReleased, contaminationProb, propWildtype, releaseTimes, maxTime, maxSize = 1000000)
{
	times = c(releaseTimes, maxTime)
	n = round(params["gamma_shape"])
	types = c("Wld", "WMel", "WAlb")
	state = createBlankStartStateVector(params, types)
	
	WMel_m = round(Wld_m*(1 - propWildtype))
	WMel_f = round(Wld_f*(1 - propWildtype))
	Wld_m = round(propWildtype*Wld_m)
	Wld_f = round(propWildtype*Wld_f)
	
	state["Wld_m"] = round(Wld_m)
	state["WMel_m"] = round(WMel_m)
	state["Wld_f_Wld"] = round(Wld_f*Wld_m/(Wld_m + WMel_m))
	state["Wld_f_WMel"] = round(Wld_f*WMel_m/(Wld_m + WMel_m))
	state["WMel_f_WMel"] = round(WMel_f*WMel_m/(Wld_m + WMel_m))
	state["WMel_f_Wld"] = round(WMel_f*Wld_m/(Wld_m + WMel_m))
	state["WAlb_m"] = round(rbinom(1, round(numReleased), (1 - contaminationProb)))
	state["WAlb_f_WAlb"] = round(numReleased) - state["WAlb_m"]
	for(i in 1:n)
	{
		thisWldName = paste0("Wld_imm_", i)
		thisWMelName = paste0("WMel_imm_", i)
		state[thisWldName] = round(propWildtype*params["K_eq"]*(params["mu_f"]*params["theta"] + params["mu_m"])/(params["gamma_rate"]*(1 + params["theta"])))
		state[thisWMelName] = round((1 - propWildtype)*params["K_eq"]*(params["mu_f"]*params["theta"] + params["mu_m"])/(params["gamma_rate"]*(1 + params["theta"])))
	}
	#print(state)
	thisRound = simulateCTMC_cpp(state, types, params, times[1], times[2], maxSize, TRUE)
	#test = simulateCTMC_cpp(state, params, types[3], 1, types)
	#return(test)
	
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
		for(i in 2:(length(times) - 1))
		{
			nextStart = state[nrow(state), ]
			numMalesReleased = round(rbinom(1, numReleased, (1 - contaminationProb)))
			numFemalesReleased = numReleased - numMalesReleased
			
			nextStart["WAlb_m"] = nextStart["WAlb_m"] + numMalesReleased
			nextStart["WAlb_f_WAlb"] = nextStart["WAlb_f_WAlb"] + numFemalesReleased
			thisRound = simulateCTMC_cpp(nextStart, types, params, times[i], times[i + 1], maxSize, TRUE)
			t = c(t, thisRound$times)
			state = rbind(state, thisRound$states)
		}
	}
	
	return(list(t = t, state = state))
}

