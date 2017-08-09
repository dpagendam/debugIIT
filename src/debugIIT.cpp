
#include <cmath>
#include <string>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List rate_m_cpp(Rcpp::NumericVector state, Rcpp::NumericVector params, Rcpp::CharacterVector maleTypePrefix)
{
	// This function works out the transition rates and the transition state change vector for males of a specified type
	int n;
	int l;
	int ind;
	std::string lastEggStateName;
	std::string maleStateName;
	std::string maleDeathRateName;
	Rcpp::NumericVector m_birth_stateChange;
	Rcpp::NumericVector m_death_stateChange;
	Rcpp::NumericVector rates = NumericVector::create(_["birth"] = 0.0, _["death"] = 0.0); 
	Rcpp::NumericMatrix stateChange(2, state.length());
	Rcpp::List RcppOutput;
	Rcpp::CharacterVector stateNames;
	
	
	n = (int) round(params["gamma_shape"]);
	l = state.length();
	
	std::stringstream ss_imm;
	ss_imm << maleTypePrefix[0] << "_imm_" << n;
	lastEggStateName = ss_imm.str();
	
	std::stringstream ss_m;
	ss_m << maleTypePrefix[0] << "_m";
	maleStateName = ss_m.str();
	
	std::stringstream ss_mu_m;
	ss_mu_m << "mu_m_" << maleTypePrefix[0];
	maleDeathRateName = ss_mu_m.str();
	
	stateNames = state.attr("names");
	
	for(int i = 0; i < l; i++)
	{
		if(stateNames[i] == lastEggStateName)
		{
			ind = i;
			break;
		}
	}
	
	double m_birth = 0.5*params["gamma_rate"]*state[ind];
	m_birth_stateChange = rep(0.0, state.length());
	m_birth_stateChange.attr("names") = state.attr("names");
	m_birth_stateChange[maleStateName] = m_birth_stateChange[maleStateName] + 1.0;
	m_birth_stateChange[ind] = m_birth_stateChange[ind] - 1.0;
	
	double m_death = params[maleDeathRateName]*state[maleStateName];
	m_death_stateChange = rep(0.0, state.length());
	m_death_stateChange.attr("names") = state.attr("names");
	m_death_stateChange[maleStateName] = m_death_stateChange[maleStateName] - 1.0;
	
	rates["birth"] = m_birth;
	rates["death"] = m_death;
	
	for(int i = 0; i < l; i++)
	{
		stateChange(0, i) = m_birth_stateChange[i];
		stateChange(1, i) = m_death_stateChange[i];
	}
	
	
	RcppOutput["rates"] = rates;
	RcppOutput["stateChange"] = stateChange;

	return(RcppOutput);
}




// [[Rcpp::export]]
Rcpp::List rate_f_cpp(Rcpp::NumericVector state, Rcpp::NumericVector params, Rcpp::CharacterVector femaleTypePrefix, Rcpp::CharacterVector mateTypeSuffix, Rcpp::CharacterVector allMateTypes)
{
	// This function works out the transition rates and the transition state change vector for females of a specified type.
	// Note that unmated females should use "Unmated" for the mateTypeSuffix
	int n;
	int l;
	int ind;
	int numMateTypes = allMateTypes.length();

	std::string lastEggStateName;
	std::string femaleStateName;
	std::string unmatedFemaleStateName;
	std::string femaleDeathRateName;
	std::string malesOfThisMatingTypeName;
	Rcpp::NumericVector f_birth_stateChange;
	Rcpp::NumericVector f_death_stateChange;
	Rcpp::NumericVector rates = NumericVector::create(_["birth"] = 0.0, _["death"] = 0.0); 
	Rcpp::NumericMatrix stateChange(2, state.length());
	Rcpp::List RcppOutput;
	Rcpp::CharacterVector stateNames;
	Rcpp::CharacterVector rateNames(2 + numMateTypes);
	
	n = (int) round(params["gamma_shape"]);
	l = state.length();
	
	std::stringstream ss_imm;
	ss_imm << femaleTypePrefix[0] << "_imm_" << n;
	lastEggStateName = ss_imm.str();
	
	std::stringstream ss_f;
	ss_f << femaleTypePrefix[0] << "_f_" << mateTypeSuffix[0];
	femaleStateName = ss_f.str();
	
	std::stringstream ss_f_unmated;
	ss_f_unmated << femaleTypePrefix[0] << "_f_Unmated";
	unmatedFemaleStateName = ss_f_unmated.str();
	
	std::stringstream ss_mu_f;
	ss_mu_f << "mu_f_" << femaleTypePrefix[0];
	femaleDeathRateName = ss_mu_f.str();

	stateNames = state.attr("names");
	
	for(int i = 0; i < l; i++)
	{
		if(stateNames[i] == lastEggStateName)
		{
			ind = i;
			break;
		}
	}
	
	double f_birth = 0.0;
	if(mateTypeSuffix[0] == "Unmated")
	{
		f_birth = 0.5*params["gamma_rate"]*state[ind];
		f_birth_stateChange = rep(0.0, state.length());
		f_birth_stateChange.attr("names") = state.attr("names");
		f_birth_stateChange[unmatedFemaleStateName] = f_birth_stateChange[unmatedFemaleStateName] + 1.0;
		f_birth_stateChange[ind] = f_birth_stateChange[ind] - 1.0;
	}
	else
	{
		std::stringstream ss_m;
		ss_m << mateTypeSuffix[0] << "_m";
		malesOfThisMatingTypeName = ss_m.str();
		double malesOfThisMateType = state[malesOfThisMatingTypeName];
	
		std::stringstream ss_c;
		ss_c << "c_" << mateTypeSuffix[0];
		if(params["DD_mating"] == 1)
		{
			f_birth = params["eta"]*params[ss_c.str()]/params["H"]*state[unmatedFemaleStateName]*state[malesOfThisMatingTypeName];
		}
		else
		{
			f_birth = params["eta"]*params[ss_c.str()]/params["H"]*state[unmatedFemaleStateName];
		}
		
		f_birth_stateChange = rep(0.0, state.length());
		f_birth_stateChange.attr("names") = state.attr("names");
		f_birth_stateChange[femaleStateName] = f_birth_stateChange[femaleStateName] + 1.0;
		f_birth_stateChange[unmatedFemaleStateName] = f_birth_stateChange[unmatedFemaleStateName] - 1.0;
	}
	
	
	

	double f_death = params[femaleDeathRateName]*state[femaleStateName];
	f_death_stateChange = rep(0.0, state.length());
	f_death_stateChange.attr("names") = state.attr("names");
	f_death_stateChange[femaleStateName] = f_death_stateChange[femaleStateName] - 1.0;

	rates["birth"] = f_birth;
	rates["death"] = f_death;
	
	for(int i = 0; i < l; i++)
	{
		stateChange(0, i) = f_birth_stateChange[i];
		stateChange(1, i) = f_death_stateChange[i];
	}


	RcppOutput["rates"] = rates;
	RcppOutput["stateChange"] = stateChange;

	return(RcppOutput);
}




// [[Rcpp::export]]
Rcpp::List rate_imm_cpp(Rcpp::NumericVector state, Rcpp::NumericVector params, Rcpp::CharacterVector immatureType, Rcpp::NumericVector immatureClassNumber, Rcpp::CharacterVector allImmatureTypes)
{
	// This function works out the transition rates and the transition state change vector for immatures (future adults) of a specified type
	int n;
	int l;
	int ind;
	
	n = (int) round(params["gamma_shape"]);
	l = state.length();
	
	Rcpp::NumericMatrix imm_stateChange(1, l);
	Rcpp::NumericVector rates(1); 
	Rcpp::List RcppOutput;
	std::string stateName;
	std::string  femaleStateNameForImmature_Wld;
	std::string  femaleStateNameForImmature_Wol;

	std::stringstream ss_imm;
	ss_imm << immatureType[0] << "_imm_" << (int) immatureClassNumber[0];
	stateName = ss_imm.str();

	std::stringstream ss_f_Wld;
	ss_f_Wld << immatureType[0] << "_f_Wld";
	femaleStateNameForImmature_Wld = ss_f_Wld.str();
	std::stringstream ss_f_Wol;
	ss_f_Wol << immatureType[0] << "_f_" << immatureType[0];
	femaleStateNameForImmature_Wol = ss_f_Wol.str();

	for(int i = 0 ; i < l; i++)
	{
		imm_stateChange(0, i) = 0.0;
	}
	

	Rcpp::CharacterVector stateNames = state.attr("names");
	
	for(int i = 0; i < l; i++)
	{
		if(stateNames[i] == stateName)
		{
			ind = i;
			break;
		}
	}
	double totalImm = 0.0;
	for(int i = 0; i < allImmatureTypes.length(); i++)
	{
		for(int j = 0; j < n; j++)
		{
			std::stringstream ss;
			ss << allImmatureTypes[i] << "_imm_" << (j + 1);
			totalImm = totalImm + state[ss.str()];
		}
	}
	if(immatureClassNumber[0] == 1)
	{

		//transition is a birth into this class
		if(immatureType[0] == "Wld")
		{
			rates = params["lambda"]*state["Wld_f_Wld"]*(params["I_max"] - totalImm)/(params["I_max"]);
		}
		else
		{
			rates = params["lambda"]*state[femaleStateNameForImmature_Wol]*(params["I_max"] - totalImm)/(params["I_max"]);
			rates = rates + params["lambda"]*state[femaleStateNameForImmature_Wld]*(params["I_max"] - totalImm)/(params["I_max"]);
		}
		imm_stateChange(0, ind) = imm_stateChange[ind] + 1.0;
	}
	else
	{	

		//transition is a transition between stages in this class
		rates = state[ind-1]*params["gamma_rate"];
		imm_stateChange(0, ind) = imm_stateChange(0, ind) + 1.0;
		imm_stateChange(0, ind - 1) = imm_stateChange(0, ind - 1) - 1.0;
	}
	
	RcppOutput["rates"] = rates;
	RcppOutput["stateChange"] = imm_stateChange;
	return(RcppOutput);
}


// [[Rcpp::export]]
Rcpp::List getRates_cpp(Rcpp::NumericVector state, Rcpp::NumericVector params, Rcpp::CharacterVector allPrefixTypes)
{
	// This function aggregates all of the transition rates and transition vectors for all the possible types of events that can occur in the population
	int numTypes = allPrefixTypes.length();
	int n;
	int l;
	int ind;
	n = (int) round(params["gamma_shape"]);
	l = state.length();

	// For each type of male, there are deaths and births
	// For each type of female that has been mated by each type of female there are deaths  and births
	// For each type of female there are unmated types
	// For each type there are numTypes immatures that can transition through the stages
	int numTransitions = numTypes*2 + numTypes*(numTypes + 1)*2 + numTypes*n;
	Rcpp::NumericMatrix allTransitions(numTransitions, l);
	Rcpp::NumericVector allRates(numTransitions); 
	Rcpp::List RcppOutput;
	Rcpp::CharacterVector thisMaleName(1);
	Rcpp::CharacterVector thisFemaleName(1);
	Rcpp::CharacterVector thisImmatureName(1);

	Rcpp::CharacterVector stateNames = state.attr("names");
	int counter = 0;
	for(int i = 0; i < l; i++)
	{
		for(int j = 0; j < allPrefixTypes.length(); j++)
		{
			std::stringstream ss_m;
			ss_m << allPrefixTypes[j] << "_m";
			thisMaleName[0] = ss_m.str();

			if(stateNames[i] == thisMaleName[0])
			{
				//It is a male of this type
				Rcpp::CharacterVector thisMaleType(1);
				thisMaleType[0] = allPrefixTypes[j];
				Rcpp::List rateList = rate_m_cpp(state, params, thisMaleType);
				Rcpp::NumericVector rate = rateList["rates"];
				Rcpp::NumericMatrix stateChange = rateList["stateChange"];

				for(int x = 0; x < rate.length(); x++)
				{
					allRates[counter + x] = rate[x];
					for(int y = 0; y < l; y++)
					{
						allTransitions(counter + x, y) = stateChange(x, y);
					}

				}
				counter = counter + rate.length();
				break;
			}
			
			std::stringstream ss_f_unmated;
			ss_f_unmated << allPrefixTypes[j] << "_f_Unmated";
			thisFemaleName[0] = ss_f_unmated.str();
			
			if(stateNames[i] == thisFemaleName[0])
			{
				//It is an unmated female of this type
				Rcpp::CharacterVector thisFemaleType(1);
				thisFemaleType[0] = allPrefixTypes[j];
				Rcpp::CharacterVector thisMateType(1);
				thisMateType[0] = "Unmated";
				Rcpp::List rateList = rate_f_cpp(state, params, thisFemaleType, thisMateType, allPrefixTypes);
				Rcpp::NumericVector rate = rateList["rates"];
				Rcpp::NumericMatrix stateChange = rateList["stateChange"];
				for(int x = 0; x < rate.length(); x++)
				{
					allRates[counter + x] = rate[x];
					for(int y = 0; y < l; y++)
					{
						allTransitions(counter + x, y) = stateChange(x, y);
					}
				}
				counter = counter + rate.length();
				break;
			}
			
			for(int k = 0; k < allPrefixTypes.length(); k++)
			{
				
				std::stringstream ss_f;
				ss_f << allPrefixTypes[j] << "_f_" << allPrefixTypes[k];
				thisFemaleName[0] = ss_f.str();

				if(stateNames[i] == thisFemaleName[0])
				{
					//It is a female of this type with an identified mating type
					Rcpp::CharacterVector thisFemaleType(1);
					thisFemaleType[0] = allPrefixTypes[j];
					Rcpp::CharacterVector thisMateType(1);
					thisMateType[0] = allPrefixTypes[k];
					Rcpp::List rateList = rate_f_cpp(state, params, thisFemaleType, thisMateType, allPrefixTypes);
					Rcpp::NumericVector rate = rateList["rates"];
					Rcpp::NumericMatrix stateChange = rateList["stateChange"];

					for(int x = 0; x < rate.length(); x++)
					{
						allRates[counter + x] = rate[x];
						for(int y = 0; y < l; y++)
						{
							allTransitions(counter + x, y) = stateChange(x, y);
						}

					}
					counter = counter + rate.length();
					break;
				}
			}
			
			for(int k = 0; k < n; k++)
			{
				std::stringstream ss_imm;
				ss_imm << allPrefixTypes[j] << "_imm_" << (k + 1);
				thisImmatureName[0] = ss_imm.str();

				if(stateNames[i] == thisImmatureName[0])
				{
					//Is is an immature of this type 
					Rcpp::CharacterVector thisImmatureType(1);
					thisImmatureType[0] = allPrefixTypes[j];
					Rcpp::NumericVector thisImmatureClassNumber(1);
					thisImmatureClassNumber[0] = (k + 1);

					Rcpp::List rateList = rate_imm_cpp(state, params, thisImmatureType, thisImmatureClassNumber, allPrefixTypes);

					Rcpp::NumericVector rate = rateList["rates"];
					Rcpp::NumericMatrix stateChange = rateList["stateChange"];

					for(int x = 0; x < rate.length(); x++)
					{
						allRates[counter + x] = rate[x];
						for(int y = 0; y < l; y++)
						{
							allTransitions(counter + x, y) = stateChange(x, y);
						}

					}
					counter = counter + rate.length();
					break;
				}
			}
		}
	}
	RcppOutput["rates"] = allRates;
	RcppOutput["transitions"] = allTransitions;
	return(RcppOutput);
}


// [[Rcpp::export]]
Rcpp::List simulateCTMC_cpp( Rcpp::NumericVector R_state, Rcpp::CharacterVector R_types, Rcpp::NumericVector R_params,  Rcpp::NumericVector R_startTime, Rcpp::NumericVector R_endTime, Rcpp::NumericVector R_maxSize, bool store = true)
{
	// This is the main function responsible for doing the IIT simulations
	// It gets all of the possible types of transitions that can occur and the rates at which these
	// are occurring and then simulates a continuous-time Markov chain based on this set of transition rates.
	Rcpp::NumericVector state = Rcpp::clone(R_state);
	Rcpp::NumericVector params = Rcpp::clone(R_params);
	Rcpp::NumericVector startTime = Rcpp::clone(R_startTime);
	Rcpp::NumericVector endTime = Rcpp::clone(R_endTime);
	Rcpp::CharacterVector types = Rcpp::clone(R_types);
	Rcpp::NumericVector maxSize = Rcpp::clone(R_maxSize);
	
	Rcpp::CharacterVector stateNames = state.attr("names");
	
	int n;
	n = (int) round(params["gamma_shape"]);
	int l;
	int numTypes = types.length();
	// For each type of male, there are deaths
	// For each type of female that has been mated by each type of female there are deaths 
	// For each type there are numTypes immatures that can transition through the stages
	int numTransitions = numTypes + numTypes*numTypes + numTypes*n;
	int counter = 0;
	
	
	
	l = state.length();
	Rcpp::NumericMatrix transitions(numTransitions, l);
	Rcpp::NumericVector rates(numTransitions); 
	Rcpp::List RcppOutput;

	Rcpp::NumericMatrix storeTimes(1, (int) maxSize[0]);
	Rcpp::NumericMatrix storeStates((int) maxSize[0], l);
	storeStates(0, _) = state;
	storeTimes(0, 0) = startTime[0];
	
	double cumulativeTime = startTime[0];

	int ind;
	Rcpp::NumericVector totalRate(1);

	

	while(cumulativeTime < endTime[0])
	{
		bool allZero = TRUE;
		for(int i = 0; i < l; i++)
		{
			if(state[i] != 0.0)
			{
				allZero = FALSE;
				break;
			}
		}
		if(allZero)
		{
			if(store == true)
			{
				RcppOutput["states"] = storeStates(Range(0, counter), Range(0, l - 1));
				RcppOutput["times"] = storeTimes(Range(0, 0), Range(0, counter));
				return Rcpp::wrap(RcppOutput);
			}
			else
			{
				RcppOutput["states"] = state;
				return Rcpp::wrap(RcppOutput);
			}
		}

		// getRates_cpp(Rcpp::NumericVector state, Rcpp::NumericVector params, Rcpp::CharacterVector allPrefixTypes)
		Rcpp::List transitionInfo = getRates_cpp(state, params, types);

		Rcpp::NumericVector rates = transitionInfo["rates"];
		Rcpp::NumericMatrix transitions = transitionInfo["transitions"];
		Rcpp::NumericVector cumulativeRates = cumsum(rates);
		totalRate[0] = sum(rates);

		Rcpp::NumericVector timeStep = rexp(1, totalRate[0]);

		cumulativeTime = cumulativeTime + timeStep[0];
		Rcpp::NumericVector u = runif(1)*totalRate[0];

		for(int i = -1; i < (cumulativeRates.length() - 1); i++)
		{
			if(i == -1)
			{
				if(u[0] > 0 & u[0] < cumulativeRates[0])
				{
					ind = 0;
					break;
				}
			}
			else
			{
				if(u[0] > cumulativeRates[i] & u[0] < cumulativeRates[i + 1])
				{
					ind = i + 1;
					break;
				}
			}
		}

		

		if(cumulativeTime < endTime[0])
		{
			for(int i = 0; i < l; i++)
			{
				state[i] = state[i] + transitions(ind, i);
			}
			if(store == true)
			{
				counter = counter + 1;
				storeTimes(0, counter) = cumulativeTime;
				storeStates(counter, _) = state;
			}
		
		}
		if(storeStates.nrow() == (counter + 1))
		{
			Rcpp::NumericMatrix storeTimes_temp(1, storeStates.nrow() + (int) maxSize[0]);
			Rcpp::NumericMatrix storeStates_temp(storeStates.nrow() + (int) maxSize[0], l);
			for(int i = 0; i < (counter + 1); i++)
			{
				storeTimes_temp(0, i) = storeTimes(0, i);
				for(int j = 0; j < l; j++)
				{
					storeTimes_temp(i, j) = storeStates(i, j);
				}
			}
			storeTimes = storeTimes_temp;
			storeStates = storeStates_temp;
		}
	}

	if(store == true)
	{
		RcppOutput["states"] = storeStates(Range(0, counter), Range(0, l - 1));
		RcppOutput["names"] = stateNames;
		RcppOutput["times"] = storeTimes(Range(0, 0), Range(0, counter));
		
		return Rcpp::wrap(RcppOutput);
	}
	else
	{
		RcppOutput["states"] = state;
		RcppOutput["names"] = stateNames;
		return(RcppOutput);
	}
}

// [[Rcpp::export]]
Rcpp::List getImmigrationRates_cpp(Rcpp::NumericMatrix R_blockImmigrationRates, Rcpp::CharacterVector R_types, Rcpp::NumericVector R_state)
{
	Rcpp::NumericMatrix blockImmigrationRates = Rcpp::clone(R_blockImmigrationRates);
	Rcpp::CharacterVector types = Rcpp::clone(R_types);
	Rcpp::NumericVector state = Rcpp::clone(R_state);
	
	Rcpp::List RcppOutput;
	int numBlocks = blockImmigrationRates.nrow();
	int numStates = state.size();
	int numTransitions = 0;
	int numTypes = types.size();
	// Immigration can occur for every adult type
	// numTypes of males, numTypes of unmated females, numTypes*numTypes mated female types
	int numAdultTypes = 2*numTypes + numTypes*numTypes;
	for(int i = 0; i < numBlocks; i++)
	{
		for(int j = 0; j < numBlocks; j++)
		{
			if(blockImmigrationRates(i,j) != 0.0)
			{
				numTransitions = numTransitions + numAdultTypes;
			}
		}
	}
	Rcpp::NumericMatrix transitions(numTransitions, numStates);
	Rcpp::NumericVector rates(numTransitions);
	Rcpp::CharacterVector stateNames = state.attr("names");
	Rcpp::CharacterVector thisType(2);
	Rcpp::CharacterVector thisType_short(2);
	rates = rates*0.0;
	transitions = transitions*0.0;
	
	int entryCounter = 0;
	std::stringstream ss_typeFrom;
	std::stringstream ss_typeTo;
	std::stringstream ss_typeFrom_short;
	std::stringstream ss_typeTo_short;
	int stateIndex1;
	int stateIndex2;
	bool done1 = FALSE;
	bool done2 = FALSE;
	for(int i = 0; i < numBlocks; i++)
	{
		//Block that individuals are emmigrating from
		for(int j = 0; j < numBlocks; j++)
		{
			//Block that individuals are immigrating to
			for(int x = 0; x < numTypes; x++)
			{
				//Do males
				ss_typeFrom << types[x] << "_m.Block_" << i + 1;
				thisType[0] = ss_typeFrom.str();
				ss_typeFrom_short << types[x] << "_m";
				thisType_short[0] = ss_typeFrom_short.str();
				ss_typeTo << types[x] << "_m.Block_" << j + 1;
				thisType[1] = ss_typeTo.str();
				ss_typeTo_short << types[x]  << "_m";
				thisType_short[1] = ss_typeTo_short.str();
				
				if(blockImmigrationRates(i, j) != 0.0 && thisType_short[0] == thisType_short[1])
				{
					done1 = FALSE;
					done2 = FALSE;
					for(int n = 0; n < stateNames.size(); n++)
					{
						if(stateNames(n) == thisType[0])
						{
							stateIndex1 = n;
							done1 = TRUE;
						}
						if(stateNames(n) == thisType[1])
						{
							stateIndex2 = n;
							done2 = TRUE;
						}
						if(done1 && done2)
						{
							break;
						}
					}
					rates(entryCounter) = state[stateIndex1]*blockImmigrationRates(i, j);
					transitions(entryCounter, stateIndex1) = -1;
					transitions(entryCounter, stateIndex2) = 1;
					entryCounter++;
				}
				ss_typeFrom.str("");
				ss_typeFrom.clear();
				ss_typeTo.str("");
				ss_typeTo.clear();
				ss_typeFrom_short.str("");
				ss_typeFrom_short.clear();
				ss_typeTo_short.str("");
				ss_typeTo_short.clear();

				//Do Unmated females
				ss_typeFrom << types[x] << "_f_Unmated.Block_" << i + 1;
				thisType[0] = ss_typeFrom.str();
				ss_typeFrom_short << types[x] << "_f_Unmated";
				thisType_short[0] = ss_typeFrom_short.str();
				ss_typeTo << types[x]  << "_f_Unmated.Block_" << j + 1;
				thisType[1] = ss_typeTo.str();
				ss_typeTo_short << types[x]  << "_f_Unmated";
				thisType_short[1] = ss_typeTo_short.str();
				
				
				if(blockImmigrationRates(i, j) != 0.0 && thisType_short[0] == thisType_short[1])
				{
					done1 = FALSE;
					done2 = FALSE;
					for(int n = 0; n < stateNames.size(); n++)
					{
						if(stateNames(n) == thisType[0])
						{
							stateIndex1 = n;
							done1 = TRUE;
						}
						if(stateNames(n) == thisType[1])
						{
							stateIndex2 = n;
							done2 = TRUE;
						}
						if(done1 && done2)
						{
							break;
						}
					}
					rates(entryCounter) = state[stateIndex1]*blockImmigrationRates(i, j);
					transitions(entryCounter, stateIndex1) = -1;
					transitions(entryCounter, stateIndex2) = 1;
					entryCounter++;
				}
				ss_typeFrom.str("");
				ss_typeFrom.clear();
				ss_typeTo.str("");
				ss_typeTo.clear();
				ss_typeFrom_short.str("");
				ss_typeFrom_short.clear();
				ss_typeTo_short.str("");
				ss_typeTo_short.clear();
				for(int y = 0; y < numTypes; y++)
				{
					//Do mated females
					ss_typeFrom << types[x] << "_f_" << types[y] << ".Block_" << i + 1;
					thisType[0] = ss_typeFrom.str();
					ss_typeFrom_short << types[x] << "_f_" << types[y];
					thisType_short[0] = ss_typeFrom_short.str();
					ss_typeTo << types[x]  << "_f_" << types[y] << ".Block_" << j + 1;
					thisType[1] = ss_typeTo.str();
					ss_typeTo_short << types[x] << "_f_" << types[y];
					thisType_short[1] = ss_typeTo_short.str();
					
					if(blockImmigrationRates(i, j) != 0.0 && thisType_short[0] == thisType_short[1])
					{
						done1 = FALSE;
						done2 = FALSE;
						for(int n = 0; n < stateNames.size(); n++)
						{
							if(stateNames(n) == thisType[0])
							{
								stateIndex1 = n;
								done1 = TRUE;
							}
							if(stateNames(n) == thisType[1])
							{
								stateIndex2 = n;
								done2 = TRUE;
							}
							if(done1 && done2)
							{
								break;
							}
						}
						rates(entryCounter) = state[stateIndex1]*blockImmigrationRates(i, j);
						transitions(entryCounter, stateIndex1) = -1;
						transitions(entryCounter, stateIndex2) = 1;
						entryCounter++;
					}
					ss_typeFrom.str("");
					ss_typeFrom.clear();
					ss_typeTo.str("");
					ss_typeTo.clear();
					ss_typeFrom_short.str("");
					ss_typeFrom_short.clear();
					ss_typeTo_short.str("");
					ss_typeTo_short.clear();
				}
			}
		}
	}
	RcppOutput["rates"] = rates;
	RcppOutput["transitions"] = transitions;
	return Rcpp::wrap(RcppOutput);
}


// [[Rcpp::export]]
Rcpp::List simulateMetapopulationCTMC_cpp( Rcpp::List R_stateList, Rcpp::NumericMatrix R_blockImmigrationRates, Rcpp::CharacterVector R_types, Rcpp::List R_paramList,  Rcpp::NumericVector R_startTime, Rcpp::NumericVector R_endTime, Rcpp::NumericVector maxSize, bool store = true)
{
	Rcpp::List RcppOutput;
	Rcpp::List stateList = Rcpp::clone(R_stateList);
	Rcpp::NumericMatrix blockImmigrationRates = Rcpp::clone(R_blockImmigrationRates);
	
	Rcpp::List paramList = Rcpp::clone(R_paramList);
	Rcpp::NumericVector startTime = Rcpp::clone(R_startTime);
	Rcpp::NumericVector endTime = Rcpp::clone(R_endTime);
	Rcpp::CharacterVector types = Rcpp::clone(R_types);
	
	Rcpp::List stateNamesList;
	
	Rcpp::NumericVector exampleState;
	Rcpp::NumericVector totalRate(1);
	
	Rcpp::CharacterVector thisType(1);
	
	int numBlocks = stateList.size();
	int stateLengths [numBlocks];
	int fromIndices [numBlocks];
	int toIndices [numBlocks];
	
	int numCols = 0;
	int numTransitions;
	int numTypes = types.size();
	int n;
	for(int i = 0; i < numBlocks; i++)
	{
		exampleState = stateList[i];
		stateLengths[i] = exampleState.size();
		if(i == 0)
		{
			fromIndices[i] = 0;
			toIndices[i] = stateLengths[i] - 1;
		}
		else
		{
			fromIndices[i] = toIndices[i - 1] + 1;
			toIndices[i] = fromIndices[i] + stateLengths[i] - 1;
		}
		numCols += exampleState.size();
		Rcpp::NumericVector thisParams = paramList[i];
		n = (int) round(thisParams["gamma_shape"]);
		numTransitions += numTypes*2 + numTypes*(numTypes + 1)*2 + numTypes*n;
	}
	Rcpp::NumericVector state(numCols);
	Rcpp::CharacterVector stateNames(numCols);
	for(int i = 0; i < numBlocks; i++)
	{
		Rcpp::NumericVector thisState = stateList[i];
		Rcpp::CharacterVector theseStateNames = thisState.attr("names");
		stateNamesList[i] = theseStateNames;
		for(int j = 0; j < theseStateNames.size(); j++)
		{
			std::stringstream ss_names;
			ss_names << theseStateNames[j] << ".Block_" << j;
			stateNames(fromIndices[i] + j) = ss_names.str();
		}
	}
	Rcpp::NumericMatrix storeStates_all(maxSize[0], numCols);
	Rcpp::NumericMatrix storeTimes_all(1, maxSize[0]);
	Rcpp::NumericMatrix transitions(numTransitions*numBlocks, numCols);
	Rcpp::NumericVector rates(numTransitions*numBlocks);
	double cumulativeTime = startTime[0];
	int fromIndex;
	int toIndex;
	int counter = 0;
	int colCounter = 0;
	//Put the starting state into the output matrices
	for(int i = 0; i < numBlocks; i++)
	{
		Rcpp::NumericVector thisState = stateList[i];
		Rcpp::NumericVector params = paramList[i];
		for(int j = 0; j < stateLengths[i]; j++)
		{
			state(colCounter) = thisState(j);
			storeStates_all(counter, colCounter) = thisState(j);
			colCounter++;
		}
	}
	storeTimes_all(1, counter) = cumulativeTime;
	counter++;
	//Start to simulate from the metapopulation model
	while(cumulativeTime < endTime[0])
	{
		bool allZero = TRUE;
		for(int i = 0; i < numCols; i++)
		{
			if(storeStates_all((counter - 1), i) != 0.0)
			{
				allZero = FALSE;
				break;
			}
		}
		if(allZero)
		{
			if(store == true)
			{
				RcppOutput["states"] = storeStates_all(Range(0, counter - 1), Range(0, numCols - 1));
				RcppOutput["times"] = storeTimes_all(Range(0, 0), Range(0, counter - 1));
				return Rcpp::wrap(RcppOutput);
			}
			else
			{
				RcppOutput["states"] = storeStates_all((counter - 1), _);
				return Rcpp::wrap(RcppOutput);
			}
		}
		rates = rates*0.0;
		transitions = transitions*0.0;
		int rateCounter;
		int rateLength;
		for(int i = 0; i < numBlocks; i++)
		{
			Rcpp::NumericVector thisParams = paramList[i];
			Rcpp::NumericVector thisState(stateLengths[i]);
			for(int r = fromIndices[i]; r <= toIndices[i]; r++)
			{
				thisState(r - fromIndices[i]) = storeStates_all(counter - 1, r);
			}
			Rprintf(" c3.6 ");
			thisState.attr("names"); = (Rcpp::CharacterVector) stateNamesList[i];
			Rcpp::List transitionInfo = getRates_cpp(thisState, thisParams, types);
			Rprintf(" c3.7 ");
			Rcpp::NumericVector theseRates = transitionInfo["rates"];
			rateCounter = 0;
			rateLength = theseRates.size();
			Rcpp::NumericMatrix theseTransitions = transitionInfo["transitions"];
			Rprintf(" c3.8 ");
			for(int r = 0; r < theseTransitions.nrow(); r++)
			{
				Rprintf(" c4 ");
				rates(rateCounter + r) = theseRates(r);
				for(int c = 0; c < theseTransitions.ncol(); c++)
				{
					transitions(rateCounter + r, fromIndex + c) = theseTransitions(r, c);
				}
			}
			rateCounter = rateCounter + rateLength;
		}
		
		//Get immigration rates
		Rcpp::List immigration = getImmigrationRates_cpp(blockImmigrationRates, types, state);
		Rcpp::NumericMatrix immigrationTransitions = immigration["transitions"];
		Rcpp::NumericVector immigrationRates = immigration["rates"];
		
		Rcpp::NumericVector allRates(rates.size() + immigrationRates.size());
		Rcpp::NumericMatrix allTransitions(rates.size() + immigrationRates.size(), numCols);
		
		for(int i = 0; i < rates.size(); i++)
		{
			allRates(i) = rates(i);
			for(int j = 0; j < numCols; j++)
			{
				allTransitions(i, j) = transitions(i, j);
			}
		}
		for(int i = 0; i < immigrationRates.size(); i++)
		{
			allRates(rates.size() + i) = immigrationRates(i);
			for(int j = 0; j < numCols; j++)
			{
				allTransitions(rates.size() + i, j) = immigrationTransitions(i, j);
			}
		}

		Rcpp::NumericVector cumulativeRates = cumsum(allRates);
		totalRate[0] = sum(allRates);
		
		Rcpp::NumericVector timeStep = rexp(1, totalRate[0]);

		cumulativeTime = cumulativeTime + timeStep[0];
		Rcpp::NumericVector u = runif(1)*totalRate[0];
		int ind;
		for(int i = -1; i < (cumulativeRates.length() - 1); i++)
		{
			if(i == -1)
			{
				if(u[0] > 0 & u[0] < cumulativeRates[0])
				{
					ind = 0;
					break;
				}
			}
			else
			{
				if(u[0] > cumulativeRates[i] & u[0] < cumulativeRates[i + 1])
				{
					ind = i + 1;
					break;
				}
			}
		}
		
		if(cumulativeTime < endTime[0])
		{
			for(int i = 0; i < numCols; i++)
			{
				state[i] = state[i] + allTransitions(ind, i);
			}
			if(store == true)
			{
				counter = counter + 1;
				storeTimes_all(0, counter) = cumulativeTime;
				storeStates_all(counter, _) = state;
			}
		
		}
		if(storeStates_all.nrow() == (counter + 1))
		{
			Rcpp::NumericMatrix storeTimes_temp(1, storeStates_all.nrow() + (int) maxSize[0]);
			Rcpp::NumericMatrix storeStates_temp(storeStates_all.nrow() + (int) maxSize[0], numCols);
			for(int i = 0; i < (counter + 1); i++)
			{
				storeTimes_temp(0, i) = storeTimes_all(0, i);
				for(int j = 0; j < numCols; j++)
				{
					storeTimes_temp(i, j) = storeStates_all(i, j);
				}
			}
			storeTimes_all = storeTimes_temp;
			storeStates_all = storeStates_temp;
		}
	}
	Rprintf(" D ");
	if(store == true)
	{
		RcppOutput["states"] = storeStates_all(Range(0, counter), Range(0, numCols - 1));
		RcppOutput["names"] = stateNames;
		RcppOutput["times"] = storeTimes_all(Range(0, 0), Range(0, counter));
		
		return Rcpp::wrap(RcppOutput);
	}
	else
	{
		RcppOutput["states"] = state;
		RcppOutput["names"] = stateNames;
		return Rcpp::wrap(RcppOutput);
	}
}

