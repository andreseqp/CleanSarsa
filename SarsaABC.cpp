 // ActCrit.cpp : Defines the entry point for the console application.
//
/*=============================================================================
ActCrit.cpp
===============================================================================

This is the main file of the project exploring a simple learning model in the
context of the cleaner wrasse mutualism. 
The model uses the actor-critic methods from reinforcement learning, 
to teach an agent to solve the market expertiment that cleaners face in 
experimental trials. In the market experiment individuals are offered two 
options of clients to clean. This options can be a visitor, a resident, or 
the absence of clients. The difference between the two types of 
clients is that visitors leave the cleaning station when they are not served, 
while residents wait; thus, are available in the next time step.There are two 
types of agent. Fully Aware Agents (FAA) estimate value for 
9 state-action pairs. In contrast, partially Aware agents 
(PAA) estimate value for 3 potential actions. 




Written by:

Andr�s E. Qui�ones
Posdoctoral researcher
Behavioural Ecology Group
Institute of Biology
University of Neuch�tel
Switzerland

Start date:
5 April 2017

=============================================================================*/

#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <omp.h>
// #include "tchar.h"
#include "../Cpp/Routines/C++/RandomNumbers/random.h"
//H for house pc, E for laptop, M for office
#include "../Cpp/json.hpp"       
// Header for reading and using JSON files see https://github.com/nlohmann/json


#define GET_VARIABLE_NAME(Variable) (#Variable)

using namespace std;
using json = nlohmann::json;

// General parameters

// Classes

enum client { resident, visitor, absence };																		
// clients can be resident, visitors, or be absent

enum learnScenario {nature, experiment, marketExperiment, ExtendedMarket};



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
	void act(client newOptions[], int &idNewOptions, double VisProbLeav, 
		double ResProbLeav, double VisReward, double ResReward, double inbr, 
		double outbr, learnScenario scenario);
		// function where the agent takes the action, gets reward, see new 
		//state and chooses future action
	void printIndData(ofstream &learnSeries, int &seed, double &outbr,
		double pV, double pR);
	// prints individual data from the learning process
	double getLearnPar(learPar parameter);
	// function to access the learning parameters
	void checkChoice();
	// Check that the choice taken is among one of the options, 
	//otherwise trigger an error
	void rebirth();																								
	// Function to reset private variables in an individual
	void getNewOptions(client newOptions[], int &idNewOptions, 
		double &VisProbLeav, double &ResProbLeav, 
		double &inbr, double &outbr, learnScenario &scenario);
	// Function to get new clients in the station, when in a natural environment
	void getExternalOptions(client newOptions[], int &idNewOptions, 
		double &inbr, double &outbr);		
	// After unattended clients leave or stay, get new clients
	void getExperimentalOptions();
	// Get new clients in the experimental setting
	void ObtainReward(double &ResReward, double &VisReward);
	// Get reward
	double softMax(double &value1, double &value2);
	// Calculate probability of taken a given action
	int mapOptionsDP(client options[], int &choice);			
	// default function that maps state action pairs to indexes in the array 
	//'values' where values are stored works for DPupdate and for 
	//state-action pair NOT for action estimation
	client cleanOptionsT[2];	// current cleaning options time = t
	client cleanOptionsT1[2];	// future cleaning options  time = t+1
	virtual void choice(int &StaAct1, int &StaAct2)=0;
	// Function to make a choice
	virtual int mapOptions(client options[], int &choice)=0;
	// function that maps state action pairs to indexes in the array 'values' 
	//where values are stored
	void getMarketExperiment();
	// Get new clients in the experimental setting of Olle's model
	void getExtenededMarket();
	// Get new clients in the experimental setting of Noa's experiment
	int numEst;
	// Number of estimates characterizing bhavioural options 9 for FIA
	int countExp;
	int getstate(bool time) {
		if (!time) return currStateAct;
		else return nextStateAct;
	}
	bool getChoice(bool time) {
		if (!time) return choiceT;
		else return choiceT1;
	}
protected:
	double values[9];																							
	// array storing the estimated values of state action pairs
	double DPbackup[9];		
	int DPid;
	int choiceT;// current choice 
	int choiceT1;// future choice
	double negRew_const;
	// Weight of the negative reward in the total reward obtained by an agent
	double alpha;// speed of learning
	double gamma;// importance of future rewards
	double tau;	// level of explorative behaviour. 
				//The higher, the less important values is when making decisions
	double neta;	
	// Weight of the negative reward in the total reward obtained by an agent
	double currentReward; // reward given by current state action pair
	double cumulReward;	// Cumulative reward
	int age;
	double negRew_curr; // size of the negative reward in the current time step
	double negReward;
	int currStateAct, nextStateAct;
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
	double &VisProbLeav, double &ResProbLeav, double &inbr, double &outbr, 
	learnScenario &scenario)
{
	if (choiceT == 0)		// Define the behaviour of the unattended client
	{
		if (cleanOptionsT[1] == resident)
		{
			if (rnd::uniform() > ResProbLeav) { 
				cleanOptionsT1[0] = cleanOptionsT[1], negRew_curr = 0; 
			}								
			// if the unttended client is a resident, it leaves with probability ResPropLeave
			else { negRew_curr = negRew_const; }
		}
		else if (cleanOptionsT[1] == visitor)
		{
			if (rnd::uniform() > VisProbLeav) { 
				cleanOptionsT1[0] = cleanOptionsT[1], negRew_curr = 0; 
			}
			// if the unttended client is a visitor, it leaves with probability VisPropLeave
			else { negRew_curr = negRew_const; }
		}
		else { negRew_curr = 0; }
	}
	else
	{
		if (cleanOptionsT[0] == resident)
		{
			if (rnd::uniform() > ResProbLeav) { 
				cleanOptionsT1[0] = cleanOptionsT[0], negRew_curr = 0;
			}	
			// if the unattended client is a resident, it leaves with probability ResPropLeave
			else { negRew_curr = negRew_const; }
		}
		else if (cleanOptionsT[0] == visitor)
		{
			if (rnd::uniform() > VisProbLeav) { 
				cleanOptionsT1[0] = cleanOptionsT[0], negRew_curr = 0; }
			// if the unattended client is a visitor, it leaves with probability VisPropLeave
			else { negRew_curr = negRew_const; }
		}
		else { negRew_const = 0; }
	}
	switch (scenario) {
	case nature: getExternalOptions(newOptions, idNewOptions, inbr, outbr);
		break;
	case experiment: getExperimentalOptions();
		break;
	case marketExperiment: getMarketExperiment();
		break;
	case ExtendedMarket: getExtenededMarket();
		break;
	default:cout << "unkown scenario!!" << endl;
		error("agent:getNewOptions",
			"unkown scenario");
		break;
	}
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
	else {
		cleanOptionsT1[0] = resident;
		cleanOptionsT1[1] = visitor;
		return;
	}
}

void agent::getMarketExperiment() {
	// Get new options in an experimental setting of Olle's models
	if (countExp==0){
		//(cleanOptionsT[0] == resident && cleanOptionsT[1] == visitor) {
		++countExp;
		return;
	}	// Every other option is a Resident-Visitor
	else if (countExp==1){
		/*((cleanOptionsT[0] == resident || cleanOptionsT[0] == visitor)
		&& cleanOptionsT[1] == absence) {*/
		cleanOptionsT1[0] = absence;
		cleanOptionsT1[1] = absence;
		++countExp;
		return;
	}
	else {
		cleanOptionsT1[0] = resident;
		cleanOptionsT1[1] = visitor;
		countExp = 0;
	}
}

void agent::getExtenededMarket() {
	// Get new options in the experimental setting of Noa's experiment
	if (countExp==0){
		//(cleanOptionsT[0] == resident && cleanOptionsT[1] == visitor) {
		++countExp;
		return;
	}	// Every other option is a Resident-Visitor
	else if (countExp==1){
		/*((cleanOptionsT[0] == resident || cleanOptionsT[0] == visitor)
		&& cleanOptionsT[1] == absence) {*/
		cleanOptionsT1[0] = absence;
		cleanOptionsT1[1] = absence;
		++countExp;
		return;
	}
	else {
		countExp = 0;
		double rand = rnd::uniform();
		if (rand < 0.5) {
			cleanOptionsT1[0] = resident;
			cleanOptionsT1[1] = visitor;
		}
		else if (rand<0.75) {
			cleanOptionsT1[0] = resident;
			cleanOptionsT1[1] = resident;
		}
		else {
			cleanOptionsT1[0] = visitor;
			cleanOptionsT1[1] = visitor;
		}
	}
}

void agent::act(client newOptions[], int &idNewOptions, double VisProbLeav, 
	double ResProbLeav, double VisReward, double ResReward, double inbr,
	double outbr, learnScenario scenario){
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
		 inbr, outbr, scenario);
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
	currStateAct = mapOptions(cleanOptionsT, choiceT);
	nextStateAct = mapOptions(cleanOptionsT1, choiceT1);
	values[currStateAct] += alpha*(currentReward +
		negReward*neta + gamma*values[nextStateAct] - values[currStateAct]);
}

void agent::printIndData(ofstream &learnSeries, int &seed, double &outbr, 
	double pV, double pR)
{
	learnSeries << seed << '\t' << age << '\t';
	//cout << seed << '\t' << age << '\t';
	learnSeries << alpha << '\t' << gamma << '\t';
	learnSeries << neta << '\t' << tau << '\t' << pV << '\t';
	learnSeries << pR <<  '\t' << outbr << '\t';
	learnSeries << cleanOptionsT[0] << '\t' << cleanOptionsT[1] << '\t';
	learnSeries << cleanOptionsT[choiceT] << '\t';
	//cout << cleanOptionsT[0] << '\t' << cleanOptionsT[1] << '\t' << choiceT << '\t';
	learnSeries << currentReward << '\t' << cumulReward << '\t' << negRew_curr << '\t';
	//cout << currentReward << '\t' << cumulReward << '\t';
	for (int j = 0; j < numEst; j++)
	{
		learnSeries << values[j] << '\t';
		//cout << values[j] << '\t';
	}
	learnSeries << endl;
	//cout << endl;
}

double agent::softMax(double &value1, double &value2)
{
	double prob1 = (exp(value1 / tau)) / (exp(value1 / tau) + exp(value2 / tau));
	// Calculate probability of chosing option 1
	return(prob1);
}


int agent::mapOptionsDP(client options[], int &choice){
	int state;
	if (options[0] == absence || options[1] == absence)	{
		// One of the options is empty
		if (options[0] == resident || options[1] == resident){					
		// the other one is a resident
			state = 2;                                                    // R0
		}
		else if (options[0] == visitor || options[1] == visitor){
			// the other one is a visitor
			state = 1;                                                    // V0
		}
		else { state = 5; }				                                  // 00
	}
	else if (options[0] == resident || options[1] == resident){
		// Both options have clients and one of them is a resident
		if (options[0] == visitor || options[1] == visitor){
			// the other one is a visitor
			state = 0;                                                    // RV
		}
		else { state = 3; }		                                          // RR
	}
	else { state = 4; }			                                          // VV
	return state;
}

class FAATyp1 :public agent{			// Fully Informed Agent (FIA)			
	public:
	FAATyp1(double alphaI, double gammaI, double tauI, double netaI)
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

string create_filename(std::string filename,nlohmann::json param){
	// name the file with the parameter specifications
	filename.append("ABCchain_CL");
	filename.append(itos(param["chain_length"]));
	filename.append("_sdPert");
	filename.append(douts(param["sdPert"]));
	filename.append("_seed");
	filename.append(itos(param["seed"]));
	filename.append(".txt");
	return(filename);
}

struct data_point {
	double rel_abund_clean, rel_abund_visitors, rel_abund_resid,prob_Vis_Leav;
	double market_exp_success;
};

struct model_param {
	//model_param(model_param const &obj);
	double alpha, gamma,tau, negReward;
	model_param &operator= (model_param const &rhs) {
		alpha = rhs.alpha;
		gamma = rhs.gamma;
		tau = rhs.tau;
		negReward = rhs.negReward;
		return *this;
	}
};

//model_param::model_param(model_param const &obj) {
//	alphaA = obj.alphaA;
//	alphaC = obj.alphaC;
//	gamma = obj.gamma;
//	negReward = obj.negReward;
//}

vector<data_point> read_Data(ifstream &marketData) {
	// open file
	if (!marketData.is_open()) {
		std::cerr << "error: unable to open data file\n";
		wait_for_return();
		exit(EXIT_FAILURE);
	}
	string header;
	getline(marketData,header); // skip header
	vector<data_point> data_set; 
	data_point input;
	for (;;) { // read data
		marketData >> input.rel_abund_clean;
		marketData >> input.rel_abund_visitors;
		marketData >> input.rel_abund_resid;
		marketData >> input.prob_Vis_Leav;
		marketData >> input.market_exp_success;
		data_set.emplace_back(input);
		// if end of file
		if (marketData.eof()) 	break;
	}
	return(data_set);
}

void initializeIndFile(ofstream &chainOutput,nlohmann::json param){
	std::string namedir = param["folder"];
	string IndFile = create_filename(namedir, param);
	chainOutput.open(IndFile.c_str());
	chainOutput << "iteration	" << "alpha	" << "gamma	" << "tau	" <<
		"negReward	" << "fit	" << "ratio" << endl;
}



//double calculate_fit(const std::vector< data_point>& empData,
//					 const std::vector< data_point>& simData) {
//	double fit = 0.0;
//	for(int i = 0; i < empData.size(); ++i) {
//		double d = empData[i].market_exp_success - simData[i].market_exp_success;
//		fit += d * d;
//	}			 
//	return(fit);
//}

double calculate_fit(const std::vector< data_point>& empData,
	const std::vector< data_point>& simData) {
	double sum_log_likelihood = 0.0;
	for (int i = 0; i < empData.size(); ++i) {
		sum_log_likelihood += log(simData[i].market_exp_success*empData[i].market_exp_success +
			(1 - simData[i].market_exp_success)*(1 - empData[i].market_exp_success)+
		0.0000000001);
	}
	return(sum_log_likelihood);
}

model_param perturb_parameters(model_param focal_param,json &sim_param) {
	model_param new_param;
	// also, you can throw in your own random number generator that you want,
	// I just add a number N(0, sd);
	new_param = focal_param;
	if (sim_param["pertScen"] == 0) {
		new_param.alpha += rnd::normal(0, sim_param["sdPert"]);
		clip_low(new_param.alpha, 0);
		new_param.tau +=  rnd::normal(0, sim_param["sdPert"]);
		clip_low(new_param.tau, 0);
	}
	if (sim_param["pertScen"] < 2) {
		new_param.tau += rnd::normal(0, sim_param["sdPert"]);
		clip_low(new_param.tau, 0);
		new_param.gamma +=  rnd::normal(0, sim_param["sdPert"]);
		clip_range(new_param.gamma, 0, 0.99999);
		new_param.negReward += rnd::normal(0, sim_param["sdPert"]);
		clip_low(new_param.negReward, 0);
	}
	else if (sim_param["pertScen"] == 2)
	{
		new_param.gamma += rnd::normal(0, sim_param["sdPert"]);
		clip_range(new_param.gamma, 0, 0.99999);
	}
	else
	{
		new_param.negReward + rnd::normal(0, sim_param["sdPert"]);
		clip_low(new_param.negReward, 0);
	}
	return(new_param);
}



std::vector<data_point> do_simulation(//del focal_model,
	std::vector<data_point> emp_data, model_param focal_comb,
	json sim_param) {
	client *clientSet;
	clientSet = new client[int(sim_param["totRounds"]) * 2];
	int idClientSet;
	FAATyp1 Cleaner (focal_comb.alpha, focal_comb.gamma,focal_comb.tau,
		focal_comb.negReward);
	std::vector<data_point> sim_data(emp_data.size());
	double VisPref, init;
	int countRVopt;
	// Loop through the data points
	for (int id_data_point = 0; id_data_point < emp_data.size(); ++id_data_point) {
		init = focal_comb.gamma*
			(1 - pow(1 - 
				emp_data[id_data_point].rel_abund_resid - 
				emp_data[id_data_point].rel_abund_visitors, 2)) / (1 - focal_comb.gamma);
		Cleaner.rebirth(init);
		draw(clientSet, sim_param["totRounds"], 
			emp_data[id_data_point].rel_abund_resid, 
			emp_data[id_data_point].rel_abund_visitors);
		idClientSet = 0;
		VisPref = 0, countRVopt = 0;
		// Loop through the learning rounds
		for (int trial = 0; trial < sim_param["totRounds"]; ++trial) {
			Cleaner.act(clientSet, idClientSet,	emp_data[id_data_point].prob_Vis_Leav, 
				sim_param["ResProbLeav"],sim_param["VisReward"],
				sim_param["ResReward"], sim_param["inbr"], sim_param["outbr"],
				learnScenario(sim_param["scenario"]));
			Cleaner.update();
			if (trial > int(sim_param["totRounds"]) * float(sim_param["propfullPrint"])) {
				if (Cleaner.getstate(0) ==0) {
					++countRVopt;
					if(Cleaner.cleanOptionsT[Cleaner.getChoice(0)]==visitor) ++VisPref;
				}
			}
		}
		sim_data[id_data_point].rel_abund_visitors =
			emp_data[id_data_point].rel_abund_visitors;
		sim_data[id_data_point].rel_abund_resid =
			emp_data[id_data_point].rel_abund_resid;
		sim_data[id_data_point].rel_abund_clean =
			emp_data[id_data_point].rel_abund_clean;
		sim_data[id_data_point].prob_Vis_Leav =
			emp_data[id_data_point].prob_Vis_Leav;
		if (countRVopt == 0) sim_data[id_data_point].market_exp_success = 0.5;
		else sim_data[id_data_point].market_exp_success = VisPref / countRVopt;
		Cleaner.rebirth();
	}
	delete[] clientSet;
	return(sim_data);
}


int main(int argc, char* argv[]){

	mark_time(1);

	// Hardwire parameter values:
	// Only for debugging 
	// input parameters provided by a JSON file with the following
	// structure:
	//json sim_param;
	//sim_param["totRounds"]    = 30000;
	//sim_param["ResReward"]    = 1;
	//sim_param["VisReward"]    = 1;
	//sim_param["ResProbLeav"]  = 0;
	//sim_param["scenario"]  = 0;
	//sim_param["inbr"]         = 0;
	//sim_param["outbr"]        = 0;
	//sim_param["seed"]         = 1;
	//sim_param["forRat"]       = 0.0;
	//sim_param["propfullPrint"]       = 0.7;
	//sim_param["sdPert"]       = 0.01;
	//sim_param["chain_length"]       = 1000;
	//sim_param["init"]       = {0.01, 0.01 ,0.0 ,0.0};
	//sim_param["pertScen"] = 0;
	////enum perturnScen {all,  bothFut, justGam, justNegRew};
	//sim_param["folder"]       = "I:/Projects/Clean.ActCrit/Simulations/ABCtest_/";

	ifstream marketData ("I:/Projects/Clean.ActCrit/Data/data_ABC.txt");
	


	// reading of parameters: 
	ifstream parameters(argv[1]);
	if (parameters.fail()) { cout << "JSON file failed" << endl; }
	json sim_param = nlohmann::json::parse(parameters);
	
	
	//ifstream marketData (argv[2]);



	rnd::set_seed(sim_param["seed"]);

	
	//enum model {model1, model2, model3, model4}; // Some sort of model choice
	//model focal_model = sim_param["model"]; // this would be nice
	
	vector< data_point > emp_data = read_Data(marketData); 
	// read the data

	model_param init_parameters; 
	init_parameters.alpha = sim_param["init"][0];
	init_parameters.tau = sim_param["init"][1];
	init_parameters.gamma = sim_param["init"][2];
	init_parameters.negReward = sim_param["init"][3];
	
	ofstream outfile;
	initializeIndFile(outfile,sim_param);

	// some sort of combination
	
	// we calculate the fit of the starting point, first we simulate data using 
	// the initial parameters
	// we also pass on the empirical data, to use the x and y coordinates.
 	std::vector< data_point > simulated_data = 
		do_simulation(//focal_model
			emp_data, init_parameters, sim_param);

   double fit = calculate_fit(emp_data, simulated_data); // function that calculates fit
		
	model_param focal_param = init_parameters;
	for(int r = 0; r < sim_param["chain_length"]; ++r) {  // 
		model_param new_param = perturb_parameters(focal_param,sim_param);
		simulated_data = do_simulation(//focal_model, 
			emp_data, new_param, sim_param);
		double new_fit = calculate_fit(emp_data, simulated_data); 
		double ratio = new_fit/ fit;  
		// better fit is larger, so ratio is > 1, so accept all.
		if( rnd::uniform() < ratio) {
			focal_param = new_param;
			fit = new_fit;
		}
		outfile << r << "\t";
		outfile << focal_param.alpha << "\t"
			<< focal_param.tau << "\t"
			<< focal_param.gamma << "\t"
			<< focal_param.negReward << "\t"
			<< fit << "\t";
		outfile << ratio << endl;
	}
	outfile.close();
	// done!
	return 0;
}