
#include <cmath>
#include <string>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List rate_m_cpp(Rcpp::NumericVector state, Rcpp::NumericVector params, Rcpp::CharacterVector maleTypePrefix)
{
	
	int n;
	int l;
	int ind;
	std::string lastEggStateName;
	std::string maleStateName;
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
	
	double m_death = params["mu_m"]*state[maleStateName];
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
	int n;
	int l;
	int ind;
	int numMateTypes = allMateTypes.length();

	std::string lastEggStateName;
	std::string femaleStateName;
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
	
	std::stringstream ss_m;
	ss_m << mateTypeSuffix[0] << "_m";
	malesOfThisMatingTypeName = ss_m.str();

	double allMales = 0.0;
	double malesOfThisMateType = 0.0;
	for(int i = 0; i < allMateTypes.length(); i++)
	{
		std::stringstream ss;
		ss << allMateTypes[i] << "_m";
		allMales = allMales + state(ss.str());
		if(ss.str() == malesOfThisMatingTypeName)
		{
			malesOfThisMateType = state(ss.str());
		}
	}

	stateNames = state.attr("names");
	
	for(int i = 0; i < l; i++)
	{
		if(stateNames[i] == lastEggStateName)
		{
			ind = i;
			break;
		}
	}

	double f_birth = 0.5*params["gamma_rate"]*state[ind]*malesOfThisMateType/allMales;
	f_birth_stateChange = rep(0.0, state.length());
	f_birth_stateChange.attr("names") = state.attr("names");
	f_birth_stateChange[femaleStateName] = f_birth_stateChange[femaleStateName] + 1.0;
	f_birth_stateChange[ind] = f_birth_stateChange[ind] - 1.0;

	double f_death = params["mu_f"]*state[femaleStateName];
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
	int n;
	int l;
	int ind;
	
	n = (int) round(params["gamma_shape"]);
	l = state.length();
	
	Rcpp::NumericMatrix imm_stateChange(1, l);
	Rcpp::NumericVector rates(1); 
	Rcpp::List RcppOutput;
	std::string stateName;
	std::string  maleStateNameForImmature;
	std::string  femaleStateNameForImmature_Wld;
	std::string  femaleStateNameForImmature_Wol;
	Rcpp::CharacterVector maleStateNames(allImmatureTypes.length());
	Rcpp::CharacterVector maleFriedsIndexNames(allImmatureTypes.length());

	std::stringstream ss_imm;
	ss_imm << immatureType[0] << "_imm_" << (int) immatureClassNumber[0];
	stateName = ss_imm.str();

	std::stringstream ss_m;
	ss_m << immatureType[0] << "_m";
	maleStateNameForImmature = ss_m.str();

	std::stringstream ss_f_Wld;
	ss_f_Wld << immatureType[0] << "_f_Wld";
	femaleStateNameForImmature_Wld = ss_f_Wld.str();
	std::stringstream ss_f_Wol;
	ss_f_Wol << immatureType[0] << "_f_" << immatureType[0];
	femaleStateNameForImmature_Wol = ss_f_Wol.str();

	for(int i = 0; i < maleStateNames.length(); i++)
	{
		std::stringstream ss1;
		ss1 << allImmatureTypes[i] << "_m";
		maleStateNames[i] = ss1.str();
		
		std::stringstream ss2;
		ss2 << "c_" << allImmatureTypes[i];
		maleFriedsIndexNames[i] = ss2.str();
	}

	for(int i = 0 ; i < l; i++)
	{
		imm_stateChange(0, i) = 0.0;
	}
	

	Rcpp::CharacterVector stateNames = state.attr("names");
	double numerator_thisType = 0.0;
	double denominator = 0.0;
	double numerator_wild = state["Wld_m"]*params["c_Wld"];
	for(int i = 0; i < maleStateNames.length(); i++)
	{
		for(int j = 0; j < state.length(); j++)
		{
			if(stateNames[j] == maleStateNames[i])
			{
				denominator = denominator + state[j]*params[(std::string) maleFriedsIndexNames[j]];
				break;
			}
		}
		if(maleStateNames[i] == maleStateNameForImmature)
		{
			numerator_thisType = state[i]*params[(std::string) maleFriedsIndexNames[i]];
		}
	}
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
			if(denominator > 0.0)
			{
				rates = params["lambda"]*numerator_thisType/denominator*state["Wld_f_Wld"]*(params["N_max"] - totalImm)/(params["N_max"]);
			}
			else
			{
				rates = 0.0;
			}
		}
		else
		{
			if(denominator > 0.0)
			{
				rates = params["lambda"]*numerator_thisType/denominator*state[femaleStateNameForImmature_Wol]*(params["N_max"] - totalImm)/(params["N_max"]);
				rates = rates + params["lambda"]*numerator_wild/denominator*state[femaleStateNameForImmature_Wld]*(params["N_max"] - totalImm)/(params["N_max"]);
			}
			else
			{
				rates = 0.0;
			}
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
	int numTypes = allPrefixTypes.length();
	int n;
	int l;
	int ind;
	n = (int) round(params["gamma_shape"]);
	l = state.length();

	// For each type of male, there are deaths and births
	// For each type of female that has been mated by each type of female there are deaths  and births
	// For each type there are numTypes immatures that can transition through the stages
	int numTransitions = numTypes*2 + numTypes*numTypes*2 + numTypes*n;
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
			
			for(int k = 0; k < allPrefixTypes.length(); k++)
			{
				
				std::stringstream ss_f;
				ss_f << allPrefixTypes[j] << "_f_" << allPrefixTypes[k];
				thisFemaleName[0] = ss_f.str();

				if(stateNames[i] == thisFemaleName[0])
				{
					//Is is a female of this type with an identified mating type
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
Rcpp::List simulateCTMC_cpp( Rcpp::NumericVector R_state, Rcpp::CharacterVector R_types, Rcpp::NumericVector R_params,  Rcpp::NumericVector R_startTime, Rcpp::NumericVector R_endTime, Rcpp::NumericVector maxSize, bool store = true)
{
	Rcpp::NumericVector state = Rcpp::clone(R_state);
	Rcpp::NumericVector params = Rcpp::clone(R_params);
	Rcpp::NumericVector startTime = Rcpp::clone(R_startTime);
	Rcpp::NumericVector endTime = Rcpp::clone(R_endTime);
	Rcpp::CharacterVector types = Rcpp::clone(R_types);
	
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