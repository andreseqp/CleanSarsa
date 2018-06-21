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

Andrés E. Quiñones
Posdoctoral researcher
Behavioural Ecology Group
Institute of Biology
University of Neuchâtel
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
#include "tchar.h"
#include "M:\\Routines\\C++\\RandomNumbers\\random.h" 
//H for house pc, E for laptop, M for office
#include "D:\\quinonesa\\Dropbox\C++\\json.hpp"       
// Header for reading and using JSON files see https://github.com/nlohmann/json


#define GET_VARIABLE_NAME(Variable) (#Variable)

using namespace std;
using json = nlohmann::json;

// General parameters

// Classes

enum client { resident, visitor, absence };																		
// clients can be resident, visitors, or be absent
enum learPar { alphaPar, gammaPar, tauPar, netaPar };

class agent													// Learning agent
{
public:
	agent();												// basic contructor
	agent(double alphaI, double gammaI, double tauI, double netaI);
	// constructor providing values for the learning parameters
	~agent();																
	// destructor not really necessary
	void update();
	// function that updates the value of state-action pairs according to 
	//current reward and estimates of future values
	void act(client newOptions[], int &idNewOptions, double &VisProbLeav, 
		double &ResProbLeav, double &VisReward, double &ResReward, double &inbr, 
		double &outbr, double &negativeRew, bool &experiment);
		// function where the agent takes the action, gets reward, see new 
		//state and chooses future action
	void printIndData(ofstream &learnSeries, int &seed, double &outbr);		
	// prints individual data from the learning process
	double getLearnPar(learPar parameter);
	// function to access the learning parameters
	void checkChoice();
	// Check that the choice taken is among one of the options, 
	//otherwise trigger an error
	void rebirth();																								
	// Function to reset private variables in an individual
	void agent::getNewOptions(client newOptions[], int &idNewOptions, 
		double &VisProbLeav, double &ResProbLeav, double &negativeRew, 
		double &inbr, double &outbr, bool &experiment);
	// Function to get new clients in the station, when in a natural environment
	void agent::getExternalOptions(client newOptions[], int &idNewOptions, 
		double &inbr, double &outbr);		
	// After unattended clients leave or stay, get new clients
	void agent::getExperimentalOptions();
	// Get new clients in the experimental setting
	void ObtainReward(double &ResReward, double &VisReward);
	// Get reward
	double softMax(double &value1, double &value2);
	// Calculate probability of taken a given action
	void agent::DPupdate(double probRes, double probVis, double &VisProbLeav,
		double &ResProbLeav, double &outbr, double &ResReward, 
		double &VisReward, double &negReward, ofstream &DPdata, bool &experiment);
	// Obtain expected values from a Dynamic programing algorithm
	int mapOptionsDP(client options[], int &choice);			
	// default function that maps state action pairs to indexes in the array 
	//'values' where values are stored works for DPupdate and for 
	//state-action pair NOT for action estimation
	void agent::forget(double forRat);
	// Forgetting function: stochastic change in the estimated values
	void agent::printDPData(ofstream &DPdata, double &oubr, int &time);
	client cleanOptionsT[2];	// current cleaning options time = t
	client cleanOptionsT1[2];	// future cleaning options  time = t+1
	virtual void choice(int &StaAct1, int &StaAct2)=0;
	// Function to make a choice
	virtual int mapOptions(client options[], int &choice)=0;
	// function that maps state action pairs to indexes in the array 'values' 
	//where values are stored
	int numEst;
	// Number of estimates characterizing bhavioural options 9 for FIA
protected:
	double values[9];																							
	// array storing the estimated values of state action pairs
	double DPbackup[9];		
	int DPid;
	int choiceT;// current choice 
	int choiceT1;// future choice
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
	numEst = 9;
	for (int i = 0; i < numEst; i++) { values[i] = 0, DPbackup[i]=0; }
	alpha = 0.01, gamma = 0.5, tau = 10;								
	// Default values
	neta = 0;
	cleanOptionsT[0] = absence, cleanOptionsT[1] = absence, choiceT = 0;
	cleanOptionsT1[0] = absence, cleanOptionsT1[1] = absence, choiceT1 = 0;
	currentReward = 0, cumulReward = 0;
	age = 0;
}

agent::agent(double alphaI, double gammaI, double tauI, double netaI)
// parameterized constructor
{
	numEst = 9;
	for (int i = 0; i < numEst; i++) { values[i] = 0, DPbackup[i]=0; }
	//values[4] = 10, values[2] = 10;
	alpha = alphaI, gamma = gammaI, tau = tauI;
	neta = netaI;
	cleanOptionsT[0] = absence, cleanOptionsT[1] = absence, choiceT = 0;
	cleanOptionsT1[0] = absence, cleanOptionsT1[1] = absence, choiceT1 = 0;
	currentReward = 0, cumulReward = 0;
	age = 0;
}

void agent::rebirth()
{
	age = 0;
	cleanOptionsT[0] = absence, cleanOptionsT[1] = absence;
	cleanOptionsT1[0] = absence, cleanOptionsT1[1] = absence;
	choiceT = 0, choiceT1 = 0;
	currentReward = 0;
	cumulReward = 0;
	for (int i = 0; i < numEst; i++) { values[i] = 0, DPbackup[i]=0; }
}

agent::~agent() {}		// Destructor

void agent::checkChoice()
{
	if (choiceT > 1 )
	{
		error("agent::act", "choice is not among the options");
	}
}

double agent::getLearnPar(learPar parameter)
{
	if (parameter == alphaPar) { return alpha; }
	else if (parameter == gammaPar) { return gamma; }
	else if (parameter == tauPar) {return tau; }
	else { return neta; }
}

void agent::ObtainReward(double &ResReward, double &VisReward)
{
	if (cleanOptionsT[choiceT] == resident) { 
		currentReward = ResReward, cumulReward += ResReward; 
	}					// Obtain reward if the choice is a resident
	else if (cleanOptionsT[choiceT] == visitor) { 
		currentReward = VisReward, cumulReward += VisReward; 
	}					// Obtain reward if the choice is a visitor
	else { currentReward = 0, cumulReward += 0; }
	// No reward if there is no client in the choice
}

void agent::getNewOptions(client newOptions[], int &idNewOptions, 
	double &VisProbLeav, double &ResProbLeav, double &negativeRew, 
	double &inbr, double &outbr, bool &experiment)
{
	if (choiceT == 0)		// Define the behaviour of the unattended client
	{
		if (cleanOptionsT[1] == resident)
		{
			if (rnd::uniform() > ResProbLeav) { 
				cleanOptionsT1[0] = cleanOptionsT[1], negReward = 0; 
			}								
			// if the unttended client is a resident, it leaves with probability ResPropLeave
			else { negReward = negativeRew; }
		}
		else if (cleanOptionsT[1] == visitor)
		{
			if (rnd::uniform() > VisProbLeav) { 
				cleanOptionsT1[0] = cleanOptionsT[1], negReward = 0; 
			}
			// if the unttended client is a visitor, it leaves with probability VisPropLeave
			else { negReward = negativeRew; }
		}
		else { negReward = 0; }
	}
	else
	{
		if (cleanOptionsT[0] == resident)
		{
			if (rnd::uniform() > ResProbLeav) { 
				cleanOptionsT1[0] = cleanOptionsT[0], negReward = 0; 
			}	
			// if the unattended client is a resident, it leaves with probability ResPropLeave
			else { negReward = negativeRew; }
		}
		else if (cleanOptionsT[0] == visitor)
		{
			if (rnd::uniform() > VisProbLeav) { 
				cleanOptionsT1[0] = cleanOptionsT[0], negReward = 0; }
			// if the unattended client is a visitor, it leaves with probability VisPropLeave
			else { negReward = negativeRew; }
		}
		else { negReward = 0; }
	}
	if (experiment) { getExperimentalOptions(); }
	else { getExternalOptions(newOptions, idNewOptions, inbr, outbr); }
}

void agent::getExternalOptions(client newOptions[], int &idNewOptions, 
	double &inbr, double &outbr)
{
	if (cleanOptionsT1[0] == absence){
		// If none of the clients stayed from the previous interaction
		cleanOptionsT1[0] = newOptions[idNewOptions], ++idNewOptions;
		if (cleanOptionsT1[0] == absence){
			// If the first draw does not yield a client
			cleanOptionsT1[1] = newOptions[idNewOptions], ++idNewOptions;
			return;
		}
	}
	if (cleanOptionsT1[0] != absence){
		// Fill the second option depending on the first option
		double probs[3] = { (1 - inbr)*(1 - outbr) + inbr*outbr , 0, 0 };	
		// Define probabilities depending on parameters
		probs[1] = probs[0] + inbr*(1 - outbr);
		// First prob is of a random option	
		probs[2] = probs[1] + outbr*(1 - inbr);	
		// Second and third homophily, and heterophily respectively
		if (probs[2] != 1) error("agent:getExternalOptions", 
			"probability does not sum up to 1");
		double rand = rnd::uniform();
		if (probs[0] > rand) {	
			cleanOptionsT1[1] = newOptions[idNewOptions], ++idNewOptions;
		}						// Random
		else if (probs[1] > rand)	cleanOptionsT1[1] = cleanOptionsT1[0];												
								// homophily
		else					// heterophily
		{
			if (cleanOptionsT1[0] == resident) { cleanOptionsT1[1] = visitor; }
			else { cleanOptionsT1[1] = resident; }
		}
	}
}

void agent::getExperimentalOptions() {
// Get new options in an experimental setting
	if (cleanOptionsT[0] == resident && cleanOptionsT[1] == visitor) {
		return;	
	}	// Every other option is a Resident-Visitor
	else
	{
		cleanOptionsT1[0] = resident;
		cleanOptionsT1[1] = visitor;
		return;
	}
}

void agent::act(client newOptions[], int &idNewOptions, double &VisProbLeav, 
	double &ResProbLeav, double &VisReward, double &ResReward, double &inbr,
	double &outbr, double &negativeRew, bool &experiment){
	// taking action, obatining reward, seeing new state, choosing future action
	int StaAct1, StaAct2;
	++age;																		
	// new time step
	cleanOptionsT[0] = cleanOptionsT1[0], cleanOptionsT[1] = cleanOptionsT1[1];
	// Future state becomes current state
	choiceT = choiceT1;
	// Future action becomes current action
	checkChoice();	
	// Check that the choice is among the options
	cleanOptionsT1[0] = absence, cleanOptionsT1[1] = absence;
	// Future state is unknown
	choiceT1 = 2;
	ObtainReward(VisReward,ResReward);
	getNewOptions(newOptions, idNewOptions, VisProbLeav, ResProbLeav, 
		negativeRew, inbr, outbr,experiment);
	choiceT1 = 0;															
	// Look into the value of state action pair if option 1 is chosen
	StaAct1 = mapOptions(cleanOptionsT1, choiceT1);
	choiceT1 = 1;								
	// Look into the value of state action pair if option 1 is chosen
	StaAct2 = mapOptions(cleanOptionsT1, choiceT1);
	choice(StaAct1, StaAct2);
}

void agent::update(){																								
	// change estimated value according to current reward and estimates of future state-action pair
	int currentStAct = mapOptions(cleanOptionsT, choiceT);
	int nextStAct = mapOptions(cleanOptionsT1, choiceT1);
	values[currentStAct] += alpha*(currentReward + 
		negReward*neta + gamma*values[nextStAct] - values[currentStAct]);
}

void agent::forget(double forRat)
{
	int tmpAct = rnd::integer(numEst);
	values[tmpAct] += rnd::normal(0, forRat);
}

void agent::printIndData(ofstream &learnSeries, int &seed, double &outbr)
{
	learnSeries << seed << '\t' << age << '\t';
	//cout << seed << '\t' << age << '\t';
	learnSeries << alpha << '\t' << gamma << '\t' << tau << '\t' << neta << '\t';
	learnSeries << outbr << '\t';
	learnSeries << cleanOptionsT[0] << '\t' << cleanOptionsT[1] << '\t';
	learnSeries << cleanOptionsT[choiceT] << '\t';
	//cout << cleanOptionsT[0] << '\t' << cleanOptionsT[1] << '\t' << choiceT << '\t';
	learnSeries << currentReward << '\t' << cumulReward << '\t' << negReward << '\t';
	//cout << currentReward << '\t' << cumulReward << '\t';
	for (int j = 0; j < numEst; j++)
	{
		learnSeries << values[j] << '\t';
		//cout << values[j] << '\t';
	}
	learnSeries << endl;
	//cout << endl;
}

void agent::printDPData(ofstream &DPdata, double &outbr, int &time)	{							
	// Print the outcome of the DP estimation
	DPdata << time << '\t';
	DPdata << alpha << '\t' << gamma << '\t' << tau << '\t' << neta << '\t';
	DPdata << outbr << '\t';
	for (int j = 0; j < 9; j++)
	{
		DPdata << DPbackup[j] << '\t';
		//cout << values[j] << '\t';
	}
	DPdata << endl;
	//cout << endl;
}

double agent::softMax(double &value1, double &value2)
{
	double prob1 = (exp(value1 / tau)) / (exp(value1 / tau) + exp(value2 / tau));
	// Calculate probability of chosing option 1
	return(prob1);
}

int agent::mapOptionsDP(client options[], int &choice){

	int stateAction;
	if (options[0] == absence || options[1] == absence)	{
		// One of the options is empty
		if (options[0] == resident || options[1] == resident){					
		// the other one is a resident
			if (options[choice] == resident) { stateAction = 4; } // State = R0 , action = R		
			else { stateAction = 5; }					// State = R0 , action = 0		
		}
		else if (options[0] == visitor || options[1] == visitor){
			// the other one is a resident
			if (options[choice] == visitor) { stateAction = 2; } 
			// State = V0 , action = V
			else { stateAction = 3; }
			// State = V0 , action = 0
		}
		else { stateAction = 8; }				
		// State = 00 , action = 0		// the other one is empty too
	}
	else if (options[0] == resident || options[1] == resident){
		// Both options have clients and one of them is a resident
		if (options[0] == visitor || options[1] == visitor){
			// the other one is a visitor
			if (options[choice] == resident) { stateAction = 1; }
			// State = RV , action = R		
			else { stateAction = 0; }	
			// State = RV , action = V
		}
		else { stateAction = 7; }		// State = RR , action = R		
	}
	else { stateAction = 6; }			 // State = VV , action = V
	return stateAction;
}
void agent::DPupdate(double probRes, double probVis, double &VisProbLeav, 
	double &ResProbLeav, double &outbr, double &ResReward, double &VisReward, 
	double &negativeRew, ofstream &DPdata,bool &experiment){
// Expected value according to DP algorithm
	double transProb[9] = {0,0,0,0,0,0,0,0,0};	// Transition probabilities
	double sum; 
	double rewards[9] = {VisReward*(1-neta)+neta*ResProbLeav*negativeRew,			//0
						ResReward*(1-neta)+neta*VisProbLeav*negativeRew,			//1
						VisReward*(1-neta),											//2
						neta*VisProbLeav*negativeRew,								//3
						ResReward*(1-neta),											//4
						neta*ResProbLeav*negativeRew,								//5
						VisReward*(1 - neta) + neta*VisProbLeav*negativeRew,		//6
						ResReward*(1 - neta) + neta*ResProbLeav*negativeRew,		//7
						0};															//8
	if (experiment)     // In an experimental setting
	{
		for (int k = 0; k < 1000; k++)
		{
			for (int i = 0; i < 9; i++)
			{
				DPid = i;//DPid = mapOptionsDP(cleanOptionsT, choiceT);
				if (DPid == 0)
				{
					transProb[0] = 0;
					transProb[1] = 0;
					transProb[2] = 0;
					transProb[3] = 0;
					transProb[4] = softMax(DPbackup[4], DPbackup[5])*(1-ResProbLeav);
					transProb[5] = softMax(DPbackup[5], DPbackup[4])*(1 - ResProbLeav);
					transProb[6] = 0;
					transProb[7] = 0;
					transProb[8] = ResProbLeav;
				}
				else if (DPid == 1)
				{
					transProb[0] = 0;
					transProb[1] = 0;
					transProb[2] = softMax(DPbackup[2], DPbackup[3])*(1 - VisProbLeav);
					transProb[3] = softMax(DPbackup[3], DPbackup[2])*(1 - VisProbLeav);
					transProb[4] = 0;
					transProb[5] = 0;
					transProb[6] = 0;
					transProb[7] = 0;
					transProb[8] = VisProbLeav;
				}
				else
				{
					transProb[0] = softMax(DPbackup[0], DPbackup[1]);
					transProb[1] = softMax(DPbackup[1], DPbackup[0]);
					transProb[2] = 0;
					transProb[3] = 0;
					transProb[4] = 0;
					transProb[5] = 0;
					transProb[6] = 0;
					transProb[7] = 0;
					transProb[8] = 0;
				}
				sum = 0;
				for (int j = 0; j < 9; j++)
				{
					sum += transProb[j] * (rewards[DPid] + gamma*DPbackup[j]);
				}
				DPbackup[DPid] = sum;
			}
			printDPData(DPdata, outbr, k);
		}
	}
	else{						// In a natural setting
		for (int k = 0; k < 1000; k++){
			for (int i = 0; i < 9; i++){
				DPid = i;//DPid = mapOptionsDP(cleanOptionsT, choiceT);
				if (DPid == 0 || DPid == 5 || DPid == 7){
					transProb[0] = softMax(DPbackup[0], DPbackup[1])*
						((1 - ResProbLeav)*(probVis*(1 - outbr) + outbr) + 
						ResProbLeav * (2 * probRes*probVis*(1 - outbr) + 
						outbr*(probVis + probRes)));
					transProb[1] = softMax(DPbackup[1], DPbackup[0])*
						((1 - ResProbLeav)*(probVis*(1 - outbr) + outbr) + 
						ResProbLeav * (2 * probRes*probVis*(1 - outbr) + 
						outbr*(probVis + probRes)));
					transProb[2] = softMax(DPbackup[2], DPbackup[3])*
						ResProbLeav * probVis*(1 - probRes - probVis)*(2 - outbr);
					transProb[3] = softMax(DPbackup[3], DPbackup[2])*
						ResProbLeav * probVis*(1 - probRes - probVis)*(2 - outbr);
					transProb[4] = softMax(DPbackup[4], DPbackup[5])*
						((1 - ResProbLeav)*(1 - probRes - probVis)*(1 - outbr) 
						+ ResProbLeav * (1 - probRes - probVis)*probRes*(2 - outbr));
					transProb[5] = softMax(DPbackup[5], DPbackup[4])*
						((1 - ResProbLeav)*(1 - probRes - probVis)*(1 - outbr) 
						+ ResProbLeav * (1 - probRes - probVis)*probRes*(2 - outbr));
					transProb[6] = ResProbLeav*pow(probVis, 2)*(1-outbr);
					transProb[7] = ((1 - ResProbLeav)*probRes + 
						ResProbLeav*pow(probRes, 2))*(1-outbr);
					transProb[8] = ResProbLeav*pow((1 - probRes - probVis), 2);
				}
				else if (DPid == 1 || DPid == 3 || DPid == 6){
					transProb[0] = softMax(DPbackup[0], DPbackup[1])*
						((1 - VisProbLeav)*(probRes*(1 - outbr) + outbr) + 
						VisProbLeav * (2 * probRes*probVis*(1 - outbr) + 
						outbr*(probVis + probRes)));
					transProb[1] = softMax(DPbackup[1], DPbackup[0])*
						((1 - VisProbLeav)*(probRes*(1 - outbr) + outbr) + 
						VisProbLeav * (2 * probRes*probVis*(1 - outbr) + 
						outbr*(probVis + probRes)));
					transProb[2] = softMax(DPbackup[2], DPbackup[3])*
						((1 - VisProbLeav)*(1 - probRes - probVis)*(1 - outbr) 
						+ VisProbLeav * probVis*(1 - probRes - probVis)*(2 - outbr));
					transProb[3] = softMax(DPbackup[3], DPbackup[2])*
						((1 - VisProbLeav)*(1 - probRes - probVis)*(1 - outbr) 
						+ VisProbLeav * probVis*(1 - probRes - probVis)*(2 - outbr));
					transProb[4] = softMax(DPbackup[4], DPbackup[5])*
						VisProbLeav * (1 - probRes - probVis)*probRes*(2 - outbr);
					transProb[5] = softMax(DPbackup[5], DPbackup[4])
						*VisProbLeav * (1 - probRes - probVis)*probRes*(2 - outbr);
					transProb[6] = ((1 - VisProbLeav)*probVis + 
						VisProbLeav*pow(probVis, 2))*(1 - outbr);
					transProb[7] = VisProbLeav*pow(probRes, 2)*(1 - outbr);
					transProb[8] = VisProbLeav*pow((1 - probRes - probVis), 2);
				}
				else{
					transProb[0] = softMax(DPbackup[0], DPbackup[1]) * 
						(2 * probRes*probVis*(1 - outbr) + 
							outbr*(probVis + probRes));
					transProb[1] = softMax(DPbackup[1], DPbackup[0]) * 
						(2 * probRes*probVis*(1 - outbr) + 
							outbr*(probVis + probRes));
					transProb[2] = softMax(DPbackup[2], DPbackup[3]) * 
						probVis*(1 - probRes - probVis)*(2 - outbr);
					transProb[3] = softMax(DPbackup[3], DPbackup[2]) * 
						probVis*(1 - probRes - probVis)*(2 - outbr);
					transProb[4] = softMax(DPbackup[4], DPbackup[5]) * 
						(1 - probRes - probVis) * probRes*(2 - outbr);
					transProb[5] = softMax(DPbackup[5], DPbackup[4]) * 
						(1 - probRes - probVis) * probRes*(2 - outbr);
					transProb[6] = pow(probVis, 2)*(1 - outbr);
					transProb[7] = pow(probRes, 2)*(1 - outbr);
					transProb[8] = pow((1 - probRes - probVis), 2);
				}
				sum = 0;
				for (int j = 0; j < 9; j++){
					sum +=  transProb[j] * (rewards[DPid] + gamma*DPbackup[j]);
				}
				DPbackup[DPid] = sum;
			}
			printDPData(DPdata, outbr, k);
		}
	}
}

class FIATyp1 :public agent{			// Fully Informed Agent (FIA)			
	public:
	FIATyp1(double alphaI, double gammaI, double tauI, double netaI)
		:agent(alphaI, gammaI, tauI, netaI){
	}
	virtual void choice(int &StaAct1, int &StaAct2){
		if (cleanOptionsT1[0] != absence && cleanOptionsT1[1] != absence){
			// if there are no absences, then use desicion rule
			double tautemp = getLearnPar(tauPar);
			if (rnd::uniform() < softMax(values[StaAct1], values[StaAct2])){
				choiceT1 = 0;
			}
			else { choiceT1 = 1; }
		}
		else if (cleanOptionsT1[0] == absence){										
			// if there is an absence, then chose the other option
			choiceT1 = 1;
		}
		else{
			choiceT1 = 0;
		}
	}
	//virtual void choice(int &StaAct1, int &StaAct2)
	//{
	//	double tautemp = getLearnPar(tauPar);
	//	if (rnd::uniform() < softMax(values[StaAct1],values[StaAct2]))
	//	{ choiceT1= 0; }
	//	else { choiceT1 = 1; }
	//}	
	virtual int mapOptions(client options[], int &choice){
		return(mapOptionsDP(options, choice));
	}
};

class PIATyp1 :public agent{				// Partially Informed Agent (PIA)	
	public:
	PIATyp1(double alphaI, double gammaI, double tauI, double netaI)
	:agent(alphaI, gammaI, tauI, netaI){
		numEst = 3;
	}
	virtual void choice(int &StaAct1, int &StaAct2){
		if (cleanOptionsT1[0] != absence && cleanOptionsT1[1] != absence){												
			// if there are no absences, then use desicion rule
			if (rnd::uniform() < softMax(values[StaAct1], values[StaAct2])){
				choiceT1 = 0;
			}
			else { choiceT1 = 1; }
		}
		else if (cleanOptionsT1[0] == absence){																		
			// if there is an absence, then chose the other option
			choiceT1 = 1;
		}
		else{
			choiceT1 = 0;
		}
	}
	/*virtual void choice(int &StaAct1, int &StaAct2)
	{
		if (rnd::uniform() < softMax(values[StaAct1], values[StaAct2]))
		{
			choiceT1 = 0;
		}
		else { choiceT1 = 1; }
	}*/
	int mapOptions(client options[], int &choice){ 
		if (options[choice] == resident) { return (0); }
		else if (options[choice] == visitor) { return(1); }
		else { return(2); }
		return(options[choice]); 
	}
};

// Functions external to the agent

void draw(client trainingSet[], int rounds, double probRes, double probVis){					
	// In a natural setting draw clients according to their abundance
	double cumProbs[3] = { probRes, probRes + probVis, 1 };
	double rndNum;
	for (int i = 0; i < rounds * 2; i++){
		rndNum = rnd::uniform();
		if (rndNum < cumProbs[0]) { trainingSet[i] = resident; }
		else if (rndNum < cumProbs[1]) { trainingSet[i] = visitor; }
		else { trainingSet[i] = absence; }
	}
}

std::string itos(int j){				// turns int into string
	std::stringstream s;
	s << j;
	return s.str();
}

std::string douts(double j){			// turns double into string
	std::stringstream s;
	s << j;
	return s.str();
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
	filename.append("_pV");
	filename.append(douts(pV));
	filename.append("_pR");
	filename.append(douts(pR));
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
		folder.erase(0, 6).append("_");
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
		indOutput << "Outbr" << '\t' << "Client1" << '\t' << "Client2" << '\t';
		indOutput << "Choice" << '\t' << "Current.Reward" << '\t';
		indOutput << "Cum.Reward" << '\t' << "Neg.Reward" << '\t';

		if (learner.numEst > 3) {
			indOutput << "RV.V" << '\t' << "RV.R" << '\t' << "V0.V" << '\t';
			indOutput << "V0.0" << '\t' << "R0.R" << '\t' << "R0.0" << '\t';
			indOutput << "VV.V" << '\t' << "RR.R" << '\t' << "OO.O" << '\t';
		}
		else {
			indOutput << "Resident" << '\t' << "Visitor" << '\t';
			indOutput << "Absence" << '\t';
		}
		indOutput << endl;
	}
}


int main(int argc, _TCHAR* argv[])
{
	mark_time(1);
	ifstream input(argv[1]);
	if (input.fail()) { cout << "JSON file failed" << endl; }
	json param = nlohmann::json::parse(input);

	int const totRounds = param["totRounds"];
	double ResReward = param["ResReward"];
	double VisReward = param["VisReward"];
	double ResProbLeav = param["ResProbLeav"];
	double VisProbLeav = param["VisProbLeav"];
	double negativeRew = param["negativeRew"];
	bool experiment = param["experiment"];
	double inbr = param["inbr"];
	double outbr = param["outbr"];
	int trainingRep = param["trainingRep"];
	double alphaT = param["alphaT"];
	const int numlearn = 2;
	int printGen = param["printGen"];
	int seed = param["seed"];
	double forRat = param["forRat"];

	/*int const totRounds = 30000;
	double ResReward = 10;
	double VisReward = ResReward;
	double ResProb = 0.2;
	double VisProb = ResProb;
	double ResProbLeav = 0;
	double VisProbLeav = 1;
	double negativeRew = -10;
	bool experiment = 0;
	double inbr = 0.0;
	double outbr = 0;
	int const trainingRep = 30;//30
	double alphaT = 0.01;
	const int numlearn = 2;
	int printGen =	1;
	int seed = 9;
	double forRat = 0.0;

	double resProbRang[8] = { 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8};

	double visProbRang[8] = { 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8 };

	double outbrRang[2] = { 0.0 };

	double visProbLeavRang[5] = { 0.1, 0.25, 0.5 ,0.75 , 1 };

	double gammaT;

	double gammaRange[3] = { 0.8 };

	double tauT;

	double tauRange[1] = { 10};

	double netaT;

	double netaRange[1] = { 0.5 };*/

	rnd::set_seed(seed);

	client *clientSet;
	clientSet = new client[totRounds * 2];
	int idClientSet;

	agent *learners[numlearn];
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

							learners[0] = new FIATyp1(alphaT, *itg, *itt, *itn);
							learners[1] = new PIATyp1(alphaT, *itg, *itt, *itn);
							ofstream printTest;
							ofstream DPprint;

							for (int k = 0; k < numlearn; ++k)
							{
								initializeIndFile(printTest, *learners[k],
									param, 0, *itVisProb, *itResProb);
								for (int i = 0; i < trainingRep; i++)
								{
									draw(clientSet, totRounds, *itResProb, *itVisProb);
									idClientSet = 0;
									for (int j = 0; j < totRounds; j++)
									{
										learners[k]->act(clientSet, idClientSet,
											VisProbLeav, ResProbLeav, VisReward, ResReward,
											inbr, outbr, negativeRew, experiment);
										learners[k]->update();
										learners[k]->forget(forRat);
										if (j > totRounds*0.9)
										{
											learners[k]->printIndData(printTest, i, outbr);
										}
										else if (j%printGen == 0)
										{
											learners[k]->printIndData(printTest, i, outbr);
										}
									}
									learners[k]->rebirth();
								}
								printTest.close();
								if (k == 0) {
									initializeIndFile(DPprint, *learners[0], param, 1,
										*itVisProb, *itResProb);
									learners[k]->DPupdate(*itResProb, *itVisProb, VisProbLeav,
										ResProbLeav, outbr, ResReward, VisReward,
										negativeRew, DPprint, experiment);
									DPprint.close();
								}
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

	