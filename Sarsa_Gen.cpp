// Sarsa_Gen.cpp : Defines the entry point for the console application.
//
/*=============================================================================
Sarsa_Gen.cpp
===============================================================================

This is the main file of the project exploring a simple learning model in the
context of the cleaner wrasse mutualism. A detailed description of the model
can be found in 
"E:\Dropbox\Neuchatel\Cleaners_learning_model\\model_description_06_2_2017.pdf".
The model uses the Sarsa algorithm, from the Time Difference (TD) methods from 
reinforcement learning, to teach an agent to solve the market expertiment that 
cleaners face in experimental trials. In the market experiment individuals are 
offered two options of clients to clean. This options can be a visitor, a 
resident, or the absence of clients. The difference between the two types of 
clients is that visitors leave the cleaning station when they are not served, 
while residents wait; thus, are available in the next time step.There are two 
types of agent. Fully informed agents (FIA-StatPosTyp1_new) estimate value for 
9 state-action pairs. In contrast, partially informed agents 
(PIA-ActPosTy1_new) estimate value for 3 potential actions. 




Written by:

Andr?s E. Qui?ones
Posdoctoral researcher
Behavioural Ecology Group
Institute of Biology
University of Neuch?tel
Switzerland

Start date:
6 February 2017

Last edit date:
12 October 2017


=============================================================================*/

#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include "../Cpp/Routines/C++/RandomNumbers/random.h"
#include "../Cpp/json.hpp"
// Header for reading and using JSON files see https://github.com/nlohmann/json

#include "Client.h" 

#define GET_VARIABLE_NAME(Variable) (#Variable)


using namespace std;
using json = nlohmann::json;

// General parameters

// Classes

enum learPar { alphaPar, gammaPar, tauPar, netaPar };

class agent													// Learning agent
{
public:
	agent();												// basic contructor
	agent(double alphaI, double gammaI, double tauI, double netaI,
		int _numSti, int _numFeat, double s2RI, double initVal);
	// constructor providing values for the learning parameters
	~agent();																
	// destructor not really necessary
	void update(int attenMech=0,double maxAlpha=1);
	// function that updates the value of state-action pairs according to 
	//current reward and estimates of future values
	void act(client newOptions[], int &idNewOptions, nlohmann::json param,
		rnd::discrete_distribution visitSpProb,
		rnd::discrete_distribution residSpProb);
		// function where the agent takes the action, gets reward, see new 
		//state and chooses future action
	void printIndData(ofstream &learnSeries, int seed, json &param);
	// prints individual data from the learning process
	double getLearnPar(learPar parameter);
	// function to access the learning parameters
	void checkChoice();
	// Check that the choice taken is among one of the options, 
	//otherwise trigger an error
	void rebirth(double initVal);																								
	// Function to reset private variables in an individual
	void getNewOptions(client newOptions[], int &idNewOptions,
		nlohmann::json param, rnd::discrete_distribution visitSpProb,
		rnd::discrete_distribution residSpProb);
	// Function to get new clients in the station, when in a natural environment
	void getExternalOptions(client newOptions[], int &idNewOptions,
		nlohmann::json param,
		rnd::discrete_distribution visitSpProb,
		rnd::discrete_distribution residSpProb);
	// After unattended clinets leave or stay, get new clients
	// After unattended clients leave or stay, get new clients
	void getExperimentalOptions(nlohmann::json param,
		rnd::discrete_distribution visitSpProb,
		rnd::discrete_distribution residSpProb);
	// Get new clients in the experimental setting
	void ObtainReward();
	// Get reward
	double softMax(double &value1, double &value2);
	// Calculate probability of taken a given action
		
	// default function that maps state action pairs to indexes in the array 
	//'values' where values are stored works for DPupdate and for 
	//state-action pair NOT for action estimation
	void forget(double forRat);
	// Forgetting function: stochastic change in the estimated values
	client cleanOptionsT[2];	// current cleaning options time = t
	client cleanOptionsT1[2];	// future cleaning options  time = t+1
	virtual void value() = 0;
	// Calculates the value estimation for the current state
	virtual void choice() = 0;
	// function that maps state action pairs to indexes in the array 'values' 
	//where values are stored
	int numSti, numFeat;
	// Number of estimates characterizing bhavioural options 9 for FIA
	// Calculate new \alpha (associability) for each stimuli
	void updateAlpha(double lambda, double max,int attenMech);
protected:
	double values[10][5];																							
	// array storing the estimated values of different features
	double alphas[10][5];
	// array storing the speed of learning for different stimuli dimentions
	double beta_k[10][5];
	// array storing an estimate of uncertainty in predictions
	double h_k[10][5];
	// array storing an estimate of uncertainty in predictions introduced in 
	// Sutton 1992
	double s2_k[10][5];
	// array storing an estimate of uncertainty in predictions
	double s2R;
	// Parameter setting the level of stochasticity
	double delta; // Prediction error
	int choiceT;// current choice 
	int choiceT1;// future choice
	double valuesT[2];
	// value estimated for current state
	double valuesT1[2];
	// value estimated for future state state
	double alpha;// speed of learning
	double gamma;// importance of future rewards
	double tau;	// level of explorative behaviour. 
				//The higher, the less important values is when making decisions
	double neta;	
	// Weight of the negative reward in the total reward obtained by an agent
	double currentReward; // reward given by current state action pair
	double cumulReward;	// Cumulative reward
	int age;
	double negReward;
};

// Members of agent class

agent::agent()			// basic constructor
{
	numSti = 10;
	numFeat = 5;
	s2R = 0.25;
	for (int i = 0; i < numSti; ++i) {
		for (int j = 0; j < numFeat+1; ++j) {
			values[i][j] = 0;
			alphas[i][j] = 0.01, beta_k[i][j] = log(s2R),
				h_k[i][j] = 0, s2_k[i][j] = 0;
		}
	}
	alpha = 0.01, gamma = 0.5, tau = 10;								
	// Default values
	neta = 0,	delta = 0;
	cleanOptionsT[0] = client(), cleanOptionsT[1] = client(), choiceT = 2;
	cleanOptionsT1[0] = client(), cleanOptionsT1[1] = client(), choiceT1 = 0;
	valuesT1[0] = 0, valuesT1[1] = 0;
	currentReward = 0, cumulReward = 0;
	age = 0;
}

agent::agent(double alphaI, double gammaI, double tauI, double netaI,
	int _numSti,int _numFeat,double s2RI=0.25,double initVal = 0){
  // parameterized constructor
	numSti = _numSti, numFeat = _numFeat,s2R = s2RI;
	for (int i = 0; i < numSti; ++i) {
		for (int j = 0; j < numFeat+1; ++j) {
			values[i][j] = 0;
			alphas[i][j] = 0.01, beta_k[i][j] = log(s2R),
				h_k[i][j] = 0, s2_k[i][j] = 0;
		}
	}
	alpha = alphaI, gamma = gammaI, tau = tauI;
	neta = netaI;
	delta = 0;
	cleanOptionsT[0] = client(), cleanOptionsT[1] = client(), choiceT = 2;
	cleanOptionsT1[0] = client(), cleanOptionsT1[1] = client(), choiceT1 = 0;
	valuesT1[0] = 0, valuesT1[1] = 0;
	currentReward = 0, cumulReward = 0;
	age = 0;
}

void agent::rebirth(double initVal = 0){
	age = 0;
	cleanOptionsT[0] = client(), cleanOptionsT[1] = client(), choiceT = 2;
	cleanOptionsT1[0] = client(), cleanOptionsT1[1] = client(), choiceT1 = 0;
	valuesT1[0] = 0, valuesT1[1] = 0;
	currentReward = 0;
	cumulReward = 0;
	for (int i=0; i < numSti; ++i) {
		for (int j = 0; j < numFeat + 1; ++j)	 
			values[i][j] = initVal,	alphas[i][j] = alpha;
	}
}

agent::~agent() {}		// Destructor

void agent::checkChoice(){
	if (choiceT > 1 )	{
		error("agent::act", "choice is not among the options");
	}
}

double agent::getLearnPar(learPar parameter){
	if (parameter == alphaPar) { return alpha; }
	else if (parameter == gammaPar) { return gamma; }
	else if (parameter == tauPar) {return tau; }
	else { return neta; }
}

void agent::ObtainReward() {
	currentReward = cleanOptionsT[choiceT].reward;
	cumulReward += cleanOptionsT[choiceT].reward;
}

void agent::getNewOptions(client newOptions[], int &idNewOptions,
	nlohmann::json param, rnd::discrete_distribution visitSpProb,
	rnd::discrete_distribution residSpProb) {
	if (choiceT == 0) {
		// Define the behaviour of the unattended client
		if (cleanOptionsT[1].mytype == resident) {
			if (rnd::uniform() > param["ResProbLeav"].get<double>()) {
				// if the unttended client is a resident, 
				// it leaves with probability ResPropLeave
				cleanOptionsT1[1] = cleanOptionsT[1], negReward = 0;
			}
			else { negReward = param["negativeRew"].get<double>(); }
		}
		else if (cleanOptionsT[1].mytype == visitor) {
			if (rnd::uniform() > param["VisProbLeav"].get<double>()) {
				// if the unttended client is a visitor, 
				// it leaves with probability VisPropLeave
				cleanOptionsT1[1] = cleanOptionsT[1], negReward = 0;
			}
			else { negReward = param["negativeRew"].get<double>(); }
		}
		else { negReward = 0; }
	}
	else {
		if (cleanOptionsT[0].mytype == resident) {
			if (rnd::uniform() > param["ResProbLeav"].get<double>()) {
				// if the unattended client is a resident, 
				// it leaves with probability ResPropLeave
				cleanOptionsT1[0] = cleanOptionsT[0], negReward = 0;
			}
			else { negReward = param["negativeRew"].get<double>(); }
		}
		else if (cleanOptionsT[0].mytype == visitor) {
			if (rnd::uniform() > param["VisProbLeav"].get<double>()) {
				// if the unattended client is a visitor, 
				// it leaves with probability VisPropLeave
				cleanOptionsT1[0] = cleanOptionsT[0], negReward = 0;
			}
			else { negReward = param["negativeRew"].get<double>(); }
		}
		else { negReward = 0; }
	}
	if (param["experiment"].get<bool>()) {
		getExperimentalOptions(param,
			visitSpProb, residSpProb);
	}
	else {
		getExternalOptions(newOptions, idNewOptions, param,
			visitSpProb, residSpProb);
	}
}

void agent::getExternalOptions(client newOptions[], int &idNewOptions,
	nlohmann::json param, rnd::discrete_distribution visitSpProb,
	rnd::discrete_distribution residSpProb) {
	if (cleanOptionsT1[0].mytype + cleanOptionsT1[1].mytype > 3) {
		// If none of the clients stayed from the previous interaction
		bool randPos = rnd::bernoulli();
		cleanOptionsT1[randPos] = newOptions[idNewOptions];
		++idNewOptions;
		if (cleanOptionsT1[randPos].mytype == absence) {
			// If the first draw does not yield a client
			cleanOptionsT1[!randPos] = newOptions[idNewOptions], ++idNewOptions;
			return;
		}
	}

	if (cleanOptionsT1[0].mytype + cleanOptionsT1[1].mytype < 4) {
		// There is a client in the stations
		// Fill the alternative option depending on the available one
		bool filledPos = cleanOptionsT1[1].mytype != absence;
		double probs[3] = { (1 - static_cast<double>(param["inbr"]))*
			(1 - static_cast<double>(param["outbr"])) +
			static_cast<double>(param["inbr"])*
			static_cast<double>(param["outbr"]) , 0, 0 };
		// Define probabilities depending on parameters
		probs[1] = probs[0] + static_cast<double>(param["inbr"])*
			(1 - static_cast<double>(param["outbr"]));
		// First prob is of a random option	
		probs[2] = probs[1] + static_cast<double>(param["outbr"])*
			(1 - static_cast<double>(param["inbr"]));
		// Second and third homophily, and heterophily respectively
		if (probs[2] != 1) {
			error("agent:getExternalOptions",
				"probability does not sum up to 1");
		}
		double rand = rnd::uniform();
		if (probs[0] > rand) {
			cleanOptionsT1[!filledPos] = newOptions[idNewOptions];
			++idNewOptions;
		}
		else if (probs[1] > rand) {
			if (cleanOptionsT1[1].mytype == resident) {
				// homophily
				std::string chosenSp = "Sp";
				chosenSp.append(itos(residSpProb.sample() + 1));
				cleanOptionsT1[!filledPos].rebirth(resident,
					param["residents"][chosenSp]["alphas"],
					param["residents"][chosenSp]["betas"],
					param["residents"][chosenSp]["reward"], chosenSp);
			}
			else {
				std::string chosenSp = "Sp";
				chosenSp.append(itos(residSpProb.sample() + 1));
				cleanOptionsT1[!filledPos].rebirth(visitor,
					param["residents"][chosenSp]["alphas"],
					param["residents"][chosenSp]["betas"],
					param["residents"][chosenSp]["reward"], chosenSp);
			}
		}
		else {
			// heterophily
			if (cleanOptionsT1[0].mytype == resident) {
				std::string chosenSp = "Sp";
				chosenSp.append(itos(visitSpProb.sample() + 1));
				cleanOptionsT1[!filledPos].rebirth(visitor,
					param["residents"][chosenSp]["alphas"],
					param["residents"][chosenSp]["betas"],
					param["residents"][chosenSp]["reward"], chosenSp);
			}
			else {
				std::string chosenSp = "Sp";
				chosenSp.append(itos(residSpProb.sample() + 1));
				cleanOptionsT1[!filledPos].rebirth(resident,
					param["residents"][chosenSp]["alphas"],
					param["residents"][chosenSp]["betas"],
					param["residents"][chosenSp]["reward"], chosenSp);
			}
		}
	}
}

void agent::getExperimentalOptions(nlohmann::json param,
	rnd::discrete_distribution visitSpProb,
	rnd::discrete_distribution residSpProb) {
	// Get new options in an experimental setting
	if (cleanOptionsT[0].mytype == resident &&
		cleanOptionsT[1].mytype == visitor) {
		return;
	}
	// Every other option is a Resident-Visitor
	else {
		std::string chosenSp = "Sp";
		chosenSp.append(itos(residSpProb.sample() + 1));
		cleanOptionsT1[0].rebirth(resident,
			param["residents"][chosenSp]["alphas"],
			param["residents"][chosenSp]["betas"],
			param["residents"][chosenSp]["reward"], chosenSp);
		chosenSp = "Sp";
		chosenSp.append(itos(visitSpProb.sample() + 1));
		cleanOptionsT1[1].rebirth(visitor,
			param["residents"][chosenSp]["alphas"],
			param["residents"][chosenSp]["betas"],
			param["residents"][chosenSp]["reward"], chosenSp);
		return;
	}
}
void agent::act(client newOptions[], int &idNewOptions,
	nlohmann::json param,
	rnd::discrete_distribution visitSpProb,
	rnd::discrete_distribution residSpProb) {
	// taking action, obtaining reward, seeing new state, choosing future action
	++age;
	// new time step
	cleanOptionsT[0] = cleanOptionsT1[0], cleanOptionsT[1] = cleanOptionsT1[1];
	// Future state becomes current state
	choiceT = choiceT1;
	// Future action becomes current action
	valuesT[0] = valuesT1[0], valuesT[1] = valuesT1[1] ;
	checkChoice();
	// Check that the choice is among the options
	cleanOptionsT1[0].rebirth();
	cleanOptionsT1[1].rebirth();
	// Future state is unknown: only the first parameters matters for this function. 
	choiceT1 = 2;
	ObtainReward();
	getNewOptions(newOptions, idNewOptions, param, visitSpProb, residSpProb);
	value();
 	choice();
}

void agent::update(int attenMech, double maxAlpha){
	// change estimated value according to current reward and estimates of future state-action pair
	double lambda;
	lambda = currentReward + negReward * neta + gamma * valuesT1[choiceT1];
	delta = lambda - valuesT[choiceT];
	if (attenMech != 2) {
		for (int countStim = 0; countStim < numSti; ++countStim) {
			values[countStim][cleanOptionsT[choiceT].features[countStim]] +=
				alphas[countStim][cleanOptionsT[choiceT].features[countStim]]
				* delta;
		}
	}
	else{
		for (int countStim = 0; countStim < numSti; ++countStim) {
			values[countStim][cleanOptionsT[choiceT].features[countStim]] +=
				alphas[countStim][cleanOptionsT[choiceT].features[countStim]] 
				*lambda;
		}
	}
	updateAlpha(lambda,maxAlpha,attenMech);
}


double ReLU(double X) { return (X > 0) ? X : 0; }

void agent::updateAlpha(double lambda, double max=1, int attenMech = 0) {
  // implementation of mechanism of selective attention. Changes in
  // the speed of learning (\alpha)
	double deltaTemp, beta, sums2x2;
	switch (attenMech)	{
	case 0: // No mechanism for selective atention
			break;
	case 1:
		for (int countStim = 0; countStim < numSti; ++countStim) {
			deltaTemp = abs(lambda - valuesT[choiceT] +
				values[countStim][cleanOptionsT[choiceT].features[countStim]]) -
				//alphas[0] = alpha*(abs(lambda - values[1]) -
				abs(lambda - 
					values[countStim][cleanOptionsT[choiceT].features[countStim]]);
			if (deltaTemp != 0) 
				alphas[countStim]
				[cleanOptionsT[choiceT].features[countStim]] += alpha*deltaTemp;
			clip_range(alphas[countStim][cleanOptionsT[choiceT].features[countStim]]
				, 0, max);
			if (isnan(alphas[countStim][cleanOptionsT[choiceT].features[countStim]])) {
				wait_for_return();
			}
		}
		 
    // attention (associability) increases for good predictors
    // Based on @mackintosh_Theory_1975
		break;
	case 2:
		for (int countStim = 0; countStim < numSti; ++countStim) {
			alphas[countStim][cleanOptionsT[choiceT].features[countStim]] = 
				alpha*abs(lambda -
				values[countStim][cleanOptionsT[choiceT].features[countStim]]) +
				(1 - alpha)*alphas[countStim][cleanOptionsT[choiceT].features[countStim]];
			clip_range(alphas[countStim][cleanOptionsT[choiceT].features[countStim]], 0, max);
			if (alphas[countStim][cleanOptionsT[choiceT].features[countStim]] > 1) {
				wait_for_return();
			}
			if (isnan(alphas[countStim][cleanOptionsT[choiceT].features[countStim]])) {
				wait_for_return();
			}
		}
    //alphas[currState] = ;// attention increases with prediction error
    // Based on @pearce_Model_1980
    break;
    case 3:
		for (int countStim = 0; countStim < numSti; ++countStim) {
			deltaTemp = 
				alpha*values[countStim][cleanOptionsT[choiceT].features[countStim]];
			if (deltaTemp != 0) 
				alphas[countStim][cleanOptionsT[choiceT].features[countStim]] +=
				alpha*deltaTemp;
			clip_range(alphas[countStim][cleanOptionsT[choiceT].features[countStim]], 0, max);
			if (isnan(alphas[countStim][cleanOptionsT[choiceT].features[countStim]])) {
				wait_for_return();
			}
		}
    // attention (associability) increases for good predictors
    // Based on @mackintosh_Theory_1975
	case 4:
		for (int countStim = 0; countStim < numSti; ++countStim) {
			if (deltaTemp != 0) 
				alphas[countStim][cleanOptionsT[choiceT].features[countStim]] = alpha*delta +
				(1 - alpha)*alphas[countStim][cleanOptionsT[choiceT].features[countStim]];
			clip_range(alphas[countStim][cleanOptionsT[choiceT].features[countStim]], 0, max);
			if (isnan(alphas[countStim][cleanOptionsT[choiceT].features[countStim]])) {
				wait_for_return();
			}
		}
	  //alphas[currState] = ;// attention increases with total prediction error
	  // Based on @pearce_Model_1980
	  break;
	case 5:
		for (int countStim = 0; countStim < numSti; ++countStim) {
			deltaTemp = abs(lambda - valuesT[choiceT] +
				values[countStim][cleanOptionsT[choiceT].features[countStim]]) -
				//alphas[0] = alpha*(abs(lambda - values[1]) -
				abs(lambda - values[countStim][cleanOptionsT[choiceT].features[countStim]]);
			clip_range(deltaTemp, 0.05, 1);
			beta = abs(lambda -
				values[countStim][cleanOptionsT[choiceT].features[countStim]]) +
				(1 - alpha)*alphas[countStim][cleanOptionsT[choiceT].features[countStim]];
			clip_range(beta, 0.5, 1);
			if (deltaTemp != 0) alphas[countStim][cleanOptionsT[choiceT].features[countStim]] 
				= deltaTemp*beta;
			clip_range(alphas[countStim][cleanOptionsT[choiceT].features[countStim]], 0, max);
			if (isnan(alphas[countStim][cleanOptionsT[choiceT].features[countStim]])) {
				wait_for_return();
			}
		}
		// Based on 
		//Pearce, J. M., and N. J. Mackintosh. 2010. Two theories of attention: 
		//a review and a possible integration. Pages 11–39 in C. Mitchell 
		//and M. E. Le Pelley, eds. . Oxford University Press, Oxford.

		break;
	case 6:
		sums2x2 = 0;
		for (int countStim = 0; countStim < numSti; ++countStim) {
			beta_k[countStim][cleanOptionsT[choiceT].features[countStim]] +=
				alpha*delta*
				h_k[countStim][cleanOptionsT[choiceT].features[countStim]] * 1;
// Fix this
//cleanOptionsT[choiceT].features[countStim];  // equation (13) in Sutton 1992b
			s2_k[countStim][cleanOptionsT[choiceT].features[countStim]] = 
				exp(beta_k[countStim][cleanOptionsT[choiceT].features[countStim]]);     // equation (11) in Sutton 1992b
			sums2x2 += s2_k[countStim][cleanOptionsT[choiceT].features[countStim]] *
				1;
// Fix this
//* cleanOptionsT[choiceT].features[countStim]*
//cleanOptionsT[choiceT].features[countStim];
		}
		for (int countStim = 0; countStim < numSti; ++countStim) {
			alphas[countStim][cleanOptionsT[choiceT].features[countStim]] = 
				s2_k[countStim][cleanOptionsT[choiceT].features[countStim]] *1
// Fix this
//cleanOptionsT[choiceT].features[countStim] 
				/ (sums2x2 + s2R);// equation (10) in Sutton 1992b
			h_k[countStim][cleanOptionsT[choiceT].features[countStim]] =
				(h_k[countStim][cleanOptionsT[choiceT].features[countStim]] +
					alphas[countStim][cleanOptionsT[choiceT].features[countStim]] *
					delta)*
				ReLU(1 - alphas[countStim][cleanOptionsT[choiceT].features[countStim]] * 1);
// Fix this
//cleanOptionsT[choiceT].features[countStim]); // equation (15)
		}
		
		//error("argument out of range", CURRENT_FUNCTION);
		// Based on @dayan_Learning_2000 and @sutton_Gain_1992
		break;
	case 7:
		// // this is the IDBD algorithm of Sutton 1992a, as described in Sutton
		// // 1992b (there is a difference in notation, involving whether the cue x
		// // is part of the learning rate); note that starting value of bet needs to
		// // be right to prevent divergence of the algorithm (bet0 = log(1/stimd)
		// // according to Sutton 1992b)
		sums2x2 = 0;
		for (int countStim = 0; countStim < numSti; ++countStim) {
			beta_k[countStim][cleanOptionsT[choiceT].features[countStim]] +=
				alpha*delta*
				h_k[countStim][cleanOptionsT[choiceT].features[countStim]] * 1;
// Fix this
//cleanOptionsT[choiceT].features[countStim];  
			// equation (13) in Sutton 1992b
			alphas[countStim][cleanOptionsT[choiceT].features[countStim]] =
				exp(beta_k[countStim][cleanOptionsT[choiceT].features[countStim]]) * 1;
// Fix this
//cleanOptionsT[choiceT].features[countStim];  
			// eq. (17) in Sutton 1992b, cf (4) in Sutton 1992a
			if (alphas[countStim][cleanOptionsT[choiceT].features[countStim]] == INFINITY)
				cout << "Infinity value" << endl;
		}
		for (int countStim = 0; countStim < numSti; ++countStim) {
			// from eq. (20) in Sutton 1992b, cf (5) or fig 2 in Sutton 1992a
			h_k[countStim][cleanOptionsT[choiceT].features[countStim]] =
				h_k[countStim][cleanOptionsT[choiceT].features[countStim]] *
				ReLU(1 - alphas[countStim][cleanOptionsT[choiceT].features[countStim]] * 1)+
// Fix this
//cleanOptionsT[choiceT].features[countStim]) +
				alphas[countStim][cleanOptionsT[choiceT].features[countStim]] * 
				delta;
		}
		// Based on @dayan_Learning_2000 and @sutton_Gain_1992
		break;
  default:
    break;
  }
}


void agent::forget(double forRat) {
	for (int i=0; i < numSti;++i) {
		for (int j=0; j < numFeat; ++j) {
			if (j != cleanOptionsT[choiceT].features[i])
				values[i][j] -= forRat;
		}
	}
}

void agent::printIndData(ofstream &learnSeries, int seed, json &param) {
	learnSeries << seed << '\t' << age << '\t';
	learnSeries << alpha << '\t' << gamma << '\t' << tau << '\t' << neta << '\t';
	learnSeries << cleanOptionsT[0].mytype << '\t';
	for (int i = 0; i < numSti; ++i) learnSeries << cleanOptionsT[0].features[i] << '\t';
	learnSeries << cleanOptionsT[1].mytype << '\t';
	for (int i = 0; i < numSti; ++i) learnSeries << cleanOptionsT[1].features[i] << '\t';
	learnSeries << cleanOptionsT[choiceT].mytype << '\t';
	learnSeries << currentReward << '\t' << cumulReward << '\t' << negReward << '\t';
	for (int j = 0; j < numSti; j++){
		for (int k = 0; k < numFeat + 1; ++k) learnSeries << values[j][k] << '\t';
		learnSeries << alphas[j] << '\t';
	}
	learnSeries << endl;
}


double agent::softMax(double &value1, double &value2) {
	double prob1 = (exp(value1 / tau)) / (exp(value1 / tau) + exp(value2 / tau));
	// Calculate probability of chosing option 1
	return(prob1);
}


class PIATyp1 :public agent{				// Partially Informed Agent (PIA)	
	public:
	PIATyp1(double alphaI, double gammaI, double tauI, double netaI,
		int _numSti, int _numFeat, double initVal=0)
	:agent(alphaI, gammaI, tauI, netaI,_numSti,_numFeat){
		numSti = _numSti;
	}
	virtual void rebirth_a(double initVal=0) {
	    rebirth(initVal);
	}
	virtual void value() {
		valuesT1[0] = 0, valuesT1[1] = 0;
		for (int countStim = 0; countStim < numSti; ++countStim) {
			valuesT1[0] += values[countStim][cleanOptionsT1[0].features[countStim]];
			valuesT1[1] += values[countStim][cleanOptionsT1[1].features[countStim]];
		}

	}
	virtual void choice(){
		choiceT1 = rnd::bernoulli(softMax(valuesT1[1], valuesT1[0]));
		/*if (rnd::uniform() < softMax(valuesT1[0], valuesT1[1])){
				choiceT1 = 0;
			}
			else { choiceT1 = 1; }
		}*/
	}
	/*virtual void choice(int &StaAct1, int &StaAct2)
	{
		if (rnd::uniform() < softMax(values[StaAct1], values[StaAct2]))
		{
			choiceT1 = 0;
		}
		else { choiceT1 = 1; }
	}*/
};

// Functions external to the agent

rnd::discrete_distribution clientProbs(json param, string clientType) {
	int numSps = param[clientType].size();
	rnd::discrete_distribution SpProb(numSps);
	for (json::iterator itSpClient = param[clientType].begin();
		itSpClient != param[clientType].end(); ++itSpClient) {
		SpProb[distance(param[clientType].begin(), itSpClient)] =
			itSpClient->at("relAbun");
	}
	return (SpProb);
}

void draw(client trainingSet[], json param,
	rnd::discrete_distribution visitSpProb,
	rnd::discrete_distribution residSpProb) {
	double cumProbs[3] = { param["ResProb"].get<double>(),
		(param["ResProb"].get<double>() +
			param["VisProb"].get<double>()),	1 };
	double rndNum;
	json::iterator itSpsVis = param["visitors"].begin();
	json::iterator itSpsRes = param["residents"].begin();
	int trialSpChanges[10];
	int spNum = param["visitors"].size();
	for (int j = 0; j < param["visitors"].size() +1; ++j)
		trialSpChanges[j] = j*param["totRounds"].get<int>()*2 /
		param["visitors"].size();
	for (int i = 0; i < param["totRounds"].get<int>() * 2; i++) {
		int Sp;
		if (param["seqSp"]) {
			int j = 0;
			while (i >= trialSpChanges[j]) {
				++j, Sp = j;
			}
		}
		rndNum = rnd::uniform();
		if (rndNum < cumProbs[0]) {
			string chosenSp = "Sp";
			if (param["seqSp"]) chosenSp.append(itos(Sp));
			else chosenSp.append(itos(residSpProb.sample() + 1));
			trainingSet[i] = client(resident,
				param["residents"][chosenSp]["alphas"],
				param["residents"][chosenSp]["betas"],
				param["residents"][chosenSp]["reward"], chosenSp,
				int(param["numFeat"]));
		}
		else if (rndNum < cumProbs[1]) {
			string chosenSp = "Sp";
			if (param["seqSp"]) chosenSp.append(itos(Sp));
			else chosenSp.append(itos(residSpProb.sample() + 1));
			trainingSet[i] = client(visitor, 
				param["visitors"][chosenSp]["alphas"],
				param["visitors"][chosenSp]["betas"],
				param["visitors"][chosenSp]["reward"], chosenSp,
				int(param["numFeat"]));
		}
		else {
			trainingSet[i] = client();
		}
	}
}

string create_filename(std::string filename, agent &individual,
	nlohmann::json param, double pV, double pR) {
	// name the file with the parameter specifications
	filename.append("_alph");
	filename.append(douts(individual.getLearnPar(alphaPar)));
	filename.append("_gamma");
	filename.append(douts(individual.getLearnPar(gammaPar)));
	filename.append("_tau");
	filename.append(douts(individual.getLearnPar(tauPar)));
	filename.append("_neta");
	filename.append(douts(individual.getLearnPar(netaPar)));
	filename.append("_AttMech");
	filename.append(itos(param["attenMech"]));
	/*filename.append("_pR");
	filename.append(douts(pR));*/
	filename.append("_seed");
	filename.append(itos(param["seed"]));
	filename.append(".txt");
	return(filename);
}

void initializeIndFile(ofstream &indOutput, agent &learner,
	nlohmann::json param, bool DP, double pV, double pR) {
	std::string namedir = param["folder"];
	// "S:\\quinonesa\\Simulations\\Basic_sarsa\\"; //  //"M:\\prelim_results\\General\\"; // "E:\\Dropbox\\Neuchatel\\prelimResults\\Set_15\\IndTrain_equVal"
	std::string namedirDP = param["folder"];
	//"S:\\quinonesa\\Simulations\\Basic_sarsa\\"; //"M:\\prelim_results\\General\\"; // "E:\\Dropbox\\Neuchatel\\prelimResults\\Set_15\\IndTrain_equVal"
	std::string folder;
	if (DP){
		folder = "\\DP";
		folder.append("_");
	}
	else{
		folder = typeid(learner).name();
		folder.erase(0, 1);
		cout << folder << '\t' << learner.getLearnPar(alphaPar) << '\t';
		cout << learner.getLearnPar(gammaPar) << '\t';
		cout << learner.getLearnPar(tauPar) << '\t';
		cout << learner.getLearnPar(netaPar) << endl;
	}
	namedir.append(folder);
	string IndFile = create_filename(namedir, learner, param, pV, pR);
	indOutput.open(IndFile.c_str());
	if (DP){
		indOutput << "Time" << '\t' << "Alpha" << '\t' << "Gamma" << '\t';
		indOutput << "Tau" << '\t' << "Neta" << '\t' << "Outbr" << '\t';
		indOutput << "RV.V" << '\t' << "RV.R" << '\t' << "V0.V" << '\t';
		indOutput << "V0.0" << '\t' << "R0.R" << '\t' << "R0.0" << '\t';
		indOutput << "VV.V" << '\t' << "RR.R" << '\t' << "OO.O" << '\t';
		indOutput << endl;
	}
	else {
		indOutput << "Training" << '\t' << "Age" << '\t' << "Alpha" << '\t';
		indOutput << "Gamma" << '\t' << "Tau" << '\t' << "Neta" << '\t';
		indOutput <<  "Client1" << '\t';
		for (int k = 0; k < int(param["numSti"]); ++k)
			indOutput << "Stim1." + itos(k) << '\t';
		indOutput << "Client2" << '\t';
		for (int k = 0; k < int(param["numSti"]); ++k)
			indOutput << "Stim2." + itos(k) << '\t';
		indOutput << "Choice" << '\t' << "Current.Reward" << '\t';
		indOutput << "Cum.Reward" << '\t' << "Neg.Reward" << '\t';
		for (int i=0; i< int(param["numSti"]);++i){
			for (int j = 0; j < int(param["numFeat"]) + 1; ++j)
				indOutput << "Val." + itos(i) + itos(j) << "\t";
			indOutput << "alpha." + itos(i) << "\t";
		}
		indOutput << endl;
	}
}


int main(int argc, char* argv[])
{
	mark_time(1);

	// Hardwire parameter values:
	// Only for debugging 

	// input parameters provided by a JSON file with the following
	// structure:

	//json param;
	//param["totRounds"] = 1000;
	//param["ResReward"] = 1;
	//param["VisReward"] = param["resreward"];
	//param["ResProb"] =  0.3;
	//param["VisProb"] =  0.3;
	//param["ResProbLeav"] = 1;
	//param["VisProbLeav"] = 1;
	//param["negativeRew"] = -0.5;
	//param["experiment"] = false;
	//param["inbr"] = 0.0;
	//param["outbr"] = 0;
	//param["trainingRep"] = 10;//30
	//param["alphaT"] = 0.005;
	//param["numlearn"] = 1;
	//param["printGen"] = 1;
	//param["netaRange"] = { 0 };
	//param["gammaRange"] = { 0 };
	//param["tauRange"] = { 0.5 };
	//param["seed"] = 1;
	//param["forRat"] = 0.0;
	//param["numSti"] = 2;
	//param["numFeat"] = 2;
	//param["propfullPrint"] = 0.8;
	//param["attenMech"] = 2;
	//param["maxAlpha"] = 0.5;
	//param["seqSp"] = true;
	//param["folder"]       = "m:/Projects/CleanSarsa/Simulations/test_/";

	//param["visitors"]["Sp1"]["alphas"] = { 1 };
	//param["visitors"]["Sp1"]["betas"] = { 0.01 };
	//param["visitors"]["Sp1"]["relAbun"] = 1;
	//param["visitors"]["Sp1"]["reward"] = { 1, 0 };
	//param["visitors"]["Sp2"]["alphas"] = { 1, 0.01 };
	//param["visitors"]["Sp2"]["betas"] = { 0.01, 1 };
	//param["visitors"]["Sp2"]["relAbun"] = 1;
	//param["visitors"]["Sp2"]["reward"] = { 1, 0 };
	//param["residents"]["Sp1"]["alphas"] = { 0.5 , 1 };
	//param["residents"]["Sp1"]["betas"] = { 1,1};
	//param["residents"]["Sp1"]["relAbun"] = 1;
	//param["residents"]["Sp1"]["reward"] = { 2, 0 };
	//param["residents"]["Sp2"]["alphas"] = { 0.5 , 1 };
	//param["residents"]["Sp2"]["betas"] = { 1,1 };
	//param["residents"]["Sp2"]["relAbun"] = 1;
	//param["residents"]["Sp2"]["reward"] = { 2, 0 };
	
	//ifstream input("M:/Projects/CleanSarsa/Simulations/test_/parameters_1.json");

	// Read parameters
	ifstream input(argv[1]);
	if (input.fail()) { cout << "JSON file failed" << endl; }
	json param = nlohmann::json::parse(input);
	//
	//// Pass on parameters from JSON to c++
	//int const totRounds = param["totRounds"];
	//double ResReward = param["ResReward"];
	//double VisReward = param["VisReward"];
	//double ResProbLeav = param["ResProbLeav"];
	//double VisProbLeav = param["VisProbLeav"];
	//double negativeRew = param["negativeRew"];
	//bool experiment = param["experiment"];
	//double inbr = param["inbr"];
	//double outbr = param["outbr"];
	//int trainingRep = param["trainingRep"];
	//double alphaT = param["alphaT"];
	//int numlearn = param["numlearn"];
	//int printGen = param["printGen"];
	//double propfullPrint = param["propfullPrint"];
	//int seed = param["seed"];
	//double forRat = param["forRat"];

	

		
	rnd::set_seed(param["seed"].get<int>());
	rnd::discrete_distribution residSpProbs = clientProbs(param, "residents");
	rnd::discrete_distribution visitSpProbs = clientProbs(param, "visitors");

	client *clientSet;
	clientSet = new client[param["totRounds"].get<int>() * 2];
	int idClientSet;

	agent *learners[2];
	for (json::iterator itVisProb = param["VisProb"].begin();
		itVisProb != param["VisProb"].end(); ++itVisProb) {
		for (json::iterator itResProb = param["ResProb"].begin();
			itResProb != param["ResProb"].end(); ++itResProb) {
			double tmp1 = *itResProb;
			double tmp2 = *itVisProb;
			if (tmp1 + tmp2 <= 1) {
				for (json::iterator itn = param["netaRange"].begin();
					itn != param["netaRange"].end(); ++itn) {

					for (json::iterator itg = param["gammaRange"].begin();
						itg != param["gammaRange"].end(); ++itg) {

						for (json::iterator itt = param["tauRange"].begin();
							itt != param["tauRange"].end(); ++itt) {
							//learners[0] = new FIATyp1(alphaT, *itg, *itt, *itn);
							learners[0] = new PIATyp1(
								param["alphaT"], *itg, *itt, *itn,
								param["numSti"], param["numFeat"]);
							ofstream printTest;
							//ofstream DPprint;

							for (int k = 0; k < int(param["numlearn"]); ++k){
								initializeIndFile(printTest, *learners[k],
									param, 0, *itVisProb, *itResProb);
								for (int i = 0; i < int(param["trainingRep"]); i++){
									idClientSet = 0;
									draw(clientSet,param, visitSpProbs,
										residSpProbs);
									for (int j = 0; j < int(param["totRounds"]); j++){
										//cout << j << '\t' << endl;
										learners[k]->act(clientSet, idClientSet, param,
											visitSpProbs, residSpProbs);
										learners[k]->update(param["attenMech"],
											param["maxAlpha"]);
										learners[k]->forget(double(param["forRat"]));
										if (j > int(param["totRounds"])*double(param["propfullPrint"])){
											learners[k]->printIndData(
												printTest,i,param);
										}
										else if (j%int(param["printGen"]) == 0){
											learners[k]->printIndData(
												printTest,i, param);
										}
									}
									learners[k]->rebirth();
								}
								printTest.close();
								// if (k == 0) {
								// 	initializeIndFile(DPprint, *learners[0], param, 1,
								// 		*itVisProb, *itResProb);
								// 	learners[k]->DPupdate(*itResProb, *itVisProb, VisProbLeav,
								// 		ResProbLeav, outbr, ResReward, VisReward,
								// 		negativeRew, DPprint, experiment);
								// 	DPprint.close();
								// }
								delete learners[k];
							}

							//}
						}
					}
					//}
				}
			}
		}
	}

	delete[] clientSet;

	mark_time(0);

	//wait_for_return();

	return 0;
}

	