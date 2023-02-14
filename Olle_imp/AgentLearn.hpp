#ifndef AGLEARN_HPP
#define AGLEARN_HPP

#include "Agent.hpp"
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>


// using flt = double; // defined in Agent.hpp

//**************************** class CompStim *****************************

// This struct represents an instance of given type of a compound stimulus. An
// instance has a vector x of feature state values x[i] and a w_true of
// contributions to the reward values from the different feature states, such
// that w_true[i]*x[i] is the expected contribution from stimulus dimension i.
// The actual reward for an agent that selects the compound stimulus is the
// expected reward plus a normally distributed random 'error'.

// For our applications, which are things like foraging, with cleaner fish
// selecting between clients or bumblebee workers selecting between flowers to
// visit as examples, we think of a compound stimulus type as a species (e.g., a
// client species or a flower species). Typically, we use size as the first
// stimulus dimension, so that different species have different distributions of
// sizes (e.g., with different means and/or different variances), and we let the
// other stimulus dimensions be absence/presence. For the size dimension,
// expected reward value w_true[i] can also differ between species. For
// absence/presence dimensions, the most natural idea is that species are
// monomorphic, but that different species can either share or have different
// absence/presence features. Furthermore, the expected reward contribution
// w_true[i] can differ between species for the same trait.

// We examine cases with 4 stimulus dimensions (stimd == 4), and we study 4
// types of compound stimuli:

// cst1: first dimension variation is given by avx1 and sdx, and there is
// absence of features for the second, third, and fourth dimensions; we can
// think of this as a species of small clients.

// cst2: first dimension variation is given by avx2 and sdx, there is a feature
// present for the second dimension, and absence of feature for the third and
// fourth dimensions; we can think of this as a different species of bigger
// clients (assuming avx2 > avx1) than, in addition to size, has a feature in
// the second dimension. It can be of interest to assume that the w_true for
// this dimension is zero, so this feature can be used to discriminate species
// but it is irrelevant for reward.

// cst3: same distribution as cst2 for the first dimension, and the same
// presence of a feature in second dimension; in addition a feature present in
// the third dimension but absence in the fourth dimension. The feature in the
// third dimension has positive value (e.g., a large client species with a
// beak), making it relevant for reward.

// cst4: same distribution as cst2 for the first dimension, and the same
// presence of a feature in second dimension; in addition a feature absent in
// the third dimension but a feature present in the fourth dimension. The
// feature in the fourth dimension has negative value (e.g., some feature
// characterising a damselfish), making it relevant for reward.

template<int stimd>
struct CompStim {
public:
    using farr = Eigen::Matrix<flt, stimd, 1>;
    using rand_eng = std::mt19937;
    // using rand_uni = std::uniform_real_distribution<flt>;
    using rand_norm = std::normal_distribution<flt>;
    CompStim(flt a_sdR, int a_cstyp) :
        sdR{a_sdR},
        cstyp{a_cstyp},
        errR(0, sdR) 
        {
            w_true.fill(0); 
            x.fill(0); 
        }
    flt Reward(rand_eng& eng) { 
        return w_true.dot(x) + errR(eng);
    }
    flt sdR;       // standard deviation of the stochastic variation in reward
    int cstyp;     // indicates the type of compound stimulus
    farr w_true;   // 'true expected reward values' for the stimulus dimensions
    farr x;        // state values (cues) for the stimulus dimensions
    rand_norm errR;
};


//************************** struct LearnStat *****************************

// This struct stores data on a learning bout, allowing for different kinds of
// learning agents

template<int stimd>
struct LearnStat {
public:
    using farr = Eigen::Matrix<flt, stimd, 1>;
    int atyp;    // agent type
    int repl;    // replicate number
    int tstep;   // time step (round or trial) of learning
    int choice;  // choice (1 or 2)
    int cstyp1;  // type of compound stimulus 1
    int cstyp2;  // type of compound stimulus 2
    farr x1;     // feature values of compound stimulus 1
    farr x2;     // feature values of compound stimulus 2
    farr x;      // feature values of selected compound stimulus
    farr w_true; // true reward values for feature dimensions
    farr w;      // estimates of reward values for feature dimensions
    farr kapp;   // learning rates for features
    farr var1;   // variable for learning rates
    farr var2;   // variable for learning rates
    farr var3;   // variable for learning rates
    flt Rew;     // perceived reward
    flt Q1;      // estimated value compound stimulus 1
    flt Q2;      // estimated value compound stimulus 2
    flt Q;       // estimated value selected compound stimulus
    flt Q1true;  // true expected value compound stimulus 1
    flt Q2true;  // true expected value compound stimulus 2
    flt delt;    // 'prediction error' (TD error)
};


//*************************************************************************

// The classes below simulate learning using different meta-learning
// approaches. A learning agent experiences a number of bouts or trials. In
// each bout there is a choice between two compound stimuli, each
// characterised by feature states x[i] for one or more stimulus dimensions,
// with i = 1, ..., stimd. The feature states can either be quantitative, like
// the size of a compound stimulus, or can be the absence/presence of a
// particular feature, which is indicated as 0/1 for the corresponding x[i].
// There could be other types of stimulus dimensions, like the shape of a
// compound stimulus, but for simplicity we do not deal with those here.

// The time sequence for a bout is, first, that the agent computes the estimated
// reward values of each compound stimulus (based of the estimated values of
// feature states), then the agent selects one of the compound stimuli, and
// perceives a reward. The resulting prediction error is used to first update
// the learning rates (if the type of agent has dynamically determined learning
// rates), and then to update the estimated reward values w[i] of the features
// present in the selected compound stimulus. However, for the Pearce-Mackintosh
// algorithm, there is a different time sequence, namely to update the learning
// rates right after updating the w[i]. This latter differs from the convention
// in Sutton 1992a and 1992b (and more generally from Kalman filter updating).

// In each round an agent is faced with selecting between two compound
// stimuli, as follows. In the first T rounds, the two compound stimuli are
// independently and randomly selected to be one of type 1 and type 2. In the
// following round, there is instead a random selection between compound
// stimuli 1, 2, 3, and 4.


//**************************** class RWagLearn ****************************

// Rescorla-Wagner learning; constant and equal learning rates alph[i] for
// each stimulus dimension

template<int stimd>
class RWagLearn {
public:
    using Agent = RWagent<stimd>;
    using farr = Eigen::Matrix<flt, stimd, 1>;
    using stat_type = LearnStat<stimd>;
    using vs_type = std::vector<stat_type>;
    using comp_stim = CompStim<stimd>;
    using rand_eng = std::mt19937;
    using rand_uni = std::uniform_real_distribution<flt>;
    using rand_int = std::uniform_int_distribution<int>;
    using rand_norm = std::normal_distribution<flt>;
    using rand_discr = std::discrete_distribution<int>;
    RWagLearn(int a_atyp,
              int a_repl,
              int a_T,
              flt a_sdR,
              flt a_sdx,
              flt a_avx1,
              flt a_avx2,
              flt a_w_true1,
              flt a_w_true2,
              flt a_w_true3,
              flt a_w_true4,
              flt a_w0,
              flt a_alph0,
              flt a_omega) :
        atyp{a_atyp},
        repl{a_repl},
        T{a_T},
        sdR{a_sdR},
        sdx{a_sdx},
        avx1{a_avx1},
        avx2{a_avx2},
        w_true1{a_w_true1},
        w_true2{a_w_true2},
        w_true3{a_w_true3},
        w_true4{a_w_true4},
        w0{a_w0},
        alph0{a_alph0},
        omega{a_omega},
        agent(w0, alph0),
        cst1(sdR, 1),
        cst2(sdR, 2),
        cst3(sdR, 3),
        cst4(sdR, 4) { 
            stat.reserve(2*T); 
            Setup_CompStim();
        }
    const vs_type& Get_stat() const { return stat; }
    void Setup_CompStim();
    void Learn(rand_eng& eng);
    void Add_stat(int tstep, const comp_stim& cs1, const comp_stim& cs2);
    int atyp;        // agent type
    int repl;        // replicate number (of agent type)
    int T;           // number of learning rounds
    flt sdR;         // sd of random reward variation
    flt sdx;         // sd of variation in x for first dimension
    flt avx1;        // mean of first dim of x for cs type 1
    flt avx2;        // mean of first dim of x for cs type 2
    flt w_true1;     // reward per x for 1st stimulus dimension
    flt w_true2;     // reward for presence for 2nd dimension
    flt w_true3;     // reward for presence for 3rd dimension
    flt w_true4;     // reward for presence for 3rd dimension
    flt w0;          // starting estimated feature value
    flt alph0;       // starting feature learning rate
    flt omega;       // parameter giving choice prob from values
    Agent agent;     // learning agent
    comp_stim cst1;  // Type 1 compound stimulus
    comp_stim cst2;  // Type 2 compound stimulus
    comp_stim cst3;  // Type 3 compound stimulus
    comp_stim cst4;  // Type 4 compound stimulus
    vs_type stat;    // learning statistics
};

template<int stimd>
void RWagLearn<stimd>::Setup_CompStim()
{
    // set up compound stimulus types
    cst1.w_true[0] = w_true1;
    cst2.w_true[0] = w_true1;
    cst3.w_true[0] = w_true1;
    cst4.w_true[0] = w_true1;

    cst1.w_true[1] = w_true2;
    cst2.w_true[1] = w_true2;
    cst3.w_true[1] = w_true2;
    cst4.w_true[1] = w_true2;

    cst1.w_true[2] = w_true3;
    cst2.w_true[2] = w_true3;
    cst3.w_true[2] = w_true3;
    cst4.w_true[2] = w_true3;

    cst1.w_true[3] = w_true4;
    cst2.w_true[3] = w_true4;
    cst3.w_true[3] = w_true4;
    cst4.w_true[3] = w_true4;

    cst1.x[0] = avx1;
    cst2.x[0] = avx2;
    cst3.x[0] = avx2;
    cst4.x[0] = avx2;

    cst1.x[1] = 0;
    cst2.x[1] = 1;
    cst3.x[1] = 1;
    cst4.x[1] = 1;

    cst1.x[2] = 0;
    cst2.x[2] = 0;
    cst3.x[2] = 1;
    cst4.x[2] = 0;

    cst1.x[3] = 0;
    cst2.x[3] = 0;
    cst3.x[3] = 0;
    cst4.x[3] = 1;
}

template<int stimd>
void RWagLearn<stimd>::Learn(rand_eng& eng)
{
    // Simulate learning by Rescorla-Wagner agent over 2*T rounds. 
    rand_uni uni(0, 1);
    rand_norm nrmx(0, sdx); // random variation in first dimension
    rand_int int1_4(1, 4);  // uniform on {1, 2, 3, 4}

    // run through 2*T rounds of learning 
    for (int tstep = 1; tstep <= 2*T; ++tstep) {
        // set the two compound stimuli for this round (pedestrian!)
        comp_stim cs1 = cst1;
        comp_stim cs2 = cst2;
        if (tstep <= T) {
            cs1 = (uni(eng) < 0.5) ? cst1 : cst2;
            cs2 = (uni(eng) < 0.5) ? cst1 : cst2;
        } else {
            int typ = int1_4(eng);
            if (typ == 1) cs1 = cst1;
            else if (typ == 2) cs1 = cst2;
            else if (typ == 3) cs1 = cst3;
            else cs1 = cst4;
            typ = int1_4(eng);
            if (typ == 1) cs2 = cst1;
            else if (typ == 2) cs2 = cst2;
            else if (typ == 3) cs2 = cst3;
            else cs2 = cst4;
        }
        // and introduce random variation in x[0]
        cs1.x[0] += nrmx(eng);
        cs2.x[0] += nrmx(eng);
        // put the features states into the agent
        for (int i = 0; i < stimd; ++i) {
            agent.x1[i] = cs1.x[i];
            agent.x2[i] = cs2.x[i];
        }
        // compute estimated reward values
        agent.Q1 = agent.Qval(agent.x1);
        agent.Q2 = agent.Qval(agent.x2);
        // agent choice (based on estimated reward values)
        agent.choice = (uni(eng) < agent.PrCS1(omega)) ? 1 : 2;
        comp_stim& cs = (agent.choice == 1) ? cs1 : cs2;
        // get reward from selected compound stimulus
        agent.Rew = cs.Reward(eng);
        // store x and Q in agent
        agent.x = cs.x;
        agent.Q = (agent.choice == 1) ? agent.Q1 : agent.Q2;
        // prediction error
        agent.delt = agent.Rew - agent.Q;
        // there is no updating of learning rates for this agent
        // update estimated feature values
        agent.Updatew();
        // store learning statistics for this round
        Add_stat(tstep, cs1, cs2);
    }
}

// Add data from current round to the learning statistics container 
template<int stimd>
void RWagLearn<stimd>::Add_stat(int tstep, const comp_stim& cs1, const comp_stim& cs2) 
{
    stat_type st;
    st.atyp = atyp;
    st.repl = repl;
    st.tstep = tstep;
    st.choice = agent.choice;
    st.cstyp1 = cs1.cstyp;
    st.cstyp2 = cs2.cstyp;
    st.x1 = agent.x1;
    st.x2 = agent.x2;
    st.x = agent.x;
    st.w_true = (st.choice == 1) ? cs1.w_true : cs2.w_true; 
    st.w = agent.w;
    st.kapp.fill(0);
    st.var1 = agent.alph;;
    st.var2.fill(0);
    st.var3.fill(0);
    st.Rew = agent.Rew;
    st.Q1 = agent.Q1;
    st.Q2 = agent.Q2;
    st.Q = agent.Q;
    st.Q1true = cs1.w_true.dot(cs1.x);
    st.Q2true = cs2.w_true.dot(cs2.x);
    st.delt = agent.delt;
    stat.push_back(st);
}


//**************************** class PMagLearn ****************************

// Pearce-Mackintosh 2010 learning 

template<int stimd>
class PMagLearn {
public:
    using Agent = PMagent<stimd>;
    using farr = Eigen::Matrix<flt, stimd, 1>;
    using stat_type = LearnStat<stimd>;
    using vs_type = std::vector<stat_type>;
    using comp_stim = CompStim<stimd>;
    using rand_eng = std::mt19937;
    using rand_uni = std::uniform_real_distribution<flt>;
    using rand_int = std::uniform_int_distribution<int>;
    using rand_norm = std::normal_distribution<flt>;
    using rand_discr = std::discrete_distribution<int>;
    PMagLearn(int a_atyp,
              int a_repl,
              int a_T,
              flt a_sdR,
              flt a_sdx,
              flt a_avx1,
              flt a_avx2,
              flt a_w_true1,
              flt a_w_true2,
              flt a_w_true3,
              flt a_w_true4,
              flt a_w0,
              flt a_kapp0,
              flt a_alph0,
              flt a_sigm0,
              flt a_beta0,
              flt a_omega,
              flt a_gam_a,
              flt a_gam_s) :
        atyp{a_atyp},
        repl{a_repl},
        T{a_T},
        sdR{a_sdR},
        sdx{a_sdx},
        avx1{a_avx1},
        avx2{a_avx2},
        w_true1{a_w_true1},
        w_true2{a_w_true2},
        w_true3{a_w_true3},
        w_true4{a_w_true4},
        w0{a_w0},
        kapp0{a_kapp0},
        alph0{a_alph0},
        sigm0{a_sigm0},
        beta0{a_beta0},
        omega{a_omega},
        gam_a{a_gam_a},
        gam_s{a_gam_s},
        agent(w0, kapp0, alph0, sigm0, beta0),
        cst1(sdR, 1),
        cst2(sdR, 2),
        cst3(sdR, 3),
        cst4(sdR, 4) { 
            stat.reserve(T); 
            Setup_CompStim();
        }
    const vs_type& Get_stat() const { return stat; }
    void Setup_CompStim();
    void Learn(rand_eng& eng);
    void Add_stat(int tstep, const comp_stim& cs1, const comp_stim& cs2);
    int atyp;        // agent type
    int repl;        // replicate number (of agent type)
    int T;           // number of learning rounds
    flt sdR;         // sd of random reward variation
    flt sdx;         // sd of variation in x for first dimension
    flt avx1;        // mean of first dim of x for cs type 1
    flt avx2;        // mean of first dim of x for cs type 2
    flt w_true1;     // reward per x for 1st stimulus dimension
    flt w_true2;     // reward for presence for 2nd dimension
    flt w_true3;     // reward for presence for 3rd dimension
    flt w_true4;     // reward for presence for 3rd dimension
    flt w0;          // starting estimated feature value
    flt kapp0;       // starting feature learning rate
    flt alph0;       // starting value for alph learning element
    flt sigm0;       // starting value for sigm learning element
    flt beta0;       // starting value for beta learning element
    flt omega;       // parameter giving choice prob from values
    flt gam_a;       // meta-learning rate for alph element
    flt gam_s;       // meta-learning rate for sigm element
    Agent agent;     // learning agent
    comp_stim cst1;  // Type 1 compound stimulus
    comp_stim cst2;  // Type 2 compound stimulus
    comp_stim cst3;  // Type 3 compound stimulus
    comp_stim cst4;  // Type 4 compound stimulus
    vs_type stat;    // learning statistics
};

template<int stimd>
void PMagLearn<stimd>::Setup_CompStim()
{
    // set up compound stimulus types
    cst1.w_true[0] = w_true1;
    cst2.w_true[0] = w_true1;
    cst3.w_true[0] = w_true1;
    cst4.w_true[0] = w_true1;

    cst1.w_true[1] = w_true2;
    cst2.w_true[1] = w_true2;
    cst3.w_true[1] = w_true2;
    cst4.w_true[1] = w_true2;

    cst1.w_true[2] = w_true3;
    cst2.w_true[2] = w_true3;
    cst3.w_true[2] = w_true3;
    cst4.w_true[2] = w_true3;

    cst1.w_true[3] = w_true4;
    cst2.w_true[3] = w_true4;
    cst3.w_true[3] = w_true4;
    cst4.w_true[3] = w_true4;

    cst1.x[0] = avx1;
    cst2.x[0] = avx2;
    cst3.x[0] = avx2;
    cst4.x[0] = avx2;

    cst1.x[1] = 0;
    cst2.x[1] = 1;
    cst3.x[1] = 1;
    cst4.x[1] = 1;

    cst1.x[2] = 0;
    cst2.x[2] = 0;
    cst3.x[2] = 1;
    cst4.x[2] = 0;

    cst1.x[3] = 0;
    cst2.x[3] = 0;
    cst3.x[3] = 0;
    cst4.x[3] = 1;
}

template<int stimd>
void PMagLearn<stimd>::Learn(rand_eng& eng)
{
    // Simulate learning by Pearce-Mackintosh agent over 2*T rounds. 
    rand_uni uni(0, 1);
    rand_norm nrmx(0, sdx); // random variation in first dimension
    rand_int int1_4(1, 4);  // uniform on {1, 2, 3, 4}

    // run through 2*T rounds of learning 
    for (int tstep = 1; tstep <= 2*T; ++tstep) {
        // set the two compound stimuli for this round (pedestrian!)
        comp_stim cs1 = cst1;
        comp_stim cs2 = cst2;
        if (tstep <= T) {
            cs1 = (uni(eng) < 0.5) ? cst1 : cst2;
            cs2 = (uni(eng) < 0.5) ? cst1 : cst2;
        } else {
            int typ = int1_4(eng);
            if (typ == 1) cs1 = cst1;
            else if (typ == 2) cs1 = cst2;
            else if (typ == 3) cs1 = cst3;
            else cs1 = cst4;
            typ = int1_4(eng);
            if (typ == 1) cs2 = cst1;
            else if (typ == 2) cs2 = cst2;
            else if (typ == 3) cs2 = cst3;
            else cs2 = cst4;
        }
        // and introduce random variation in x[0]
        cs1.x[0] += nrmx(eng);
        cs2.x[0] += nrmx(eng);
        // put the features into the agent
        agent.x1 = cs1.x;
        agent.x2 = cs2.x;
        // compute estimated reward values
        agent.Q1 = agent.Qval(agent.x1);
        agent.Q2 = agent.Qval(agent.x2);
        // agent choice (based on estimated reward values)
        agent.choice = (uni(eng) < agent.PrCS1(omega)) ? 1 : 2;
        comp_stim& cs = (agent.choice == 1) ? cs1 : cs2;
        // get reward from selected compound stimulus
        agent.Rew = cs.Reward(eng);
        // store x and Q in agent
        agent.x = cs.x;
        agent.Q = (agent.choice == 1) ? agent.Q1 : agent.Q2;
        // prediction error
        agent.delt = agent.Rew - agent.Q;
        // update estimated feature values
        agent.Updatew();
        // update learning rates for the PM agent
        agent.PMupdate(gam_a, gam_s);
        // store learning statistics for this round
        Add_stat(tstep, cs1, cs2);
    }
}

// Add data from current round to the learning statistics container 
template<int stimd>
void PMagLearn<stimd>::Add_stat(int tstep, const comp_stim& cs1, const comp_stim& cs2) 
{
    stat_type st;
    st.atyp = atyp;
    st.repl = repl;
    st.tstep = tstep;
    st.choice = agent.choice;
    st.cstyp1 = cs1.cstyp;
    st.cstyp2 = cs2.cstyp;
    st.x1 = agent.x1;
    st.x2 = agent.x2;
    st.x = agent.x;
    st.w_true = (st.choice == 1) ? cs1.w_true : cs2.w_true; 
    st.w = agent.w;
    st.kapp = agent.kapp;
    st.var1 = agent.alph;
    st.var2 = agent.sigm;
    st.var3.fill(0);
    st.Rew = agent.Rew;
    st.Q1 = agent.Q1;
    st.Q2 = agent.Q2;
    st.Q = agent.Q;
    st.Q1true = cs1.w_true.dot(cs1.x);
    st.Q2true = cs2.w_true.dot(cs2.x);
    st.delt = agent.delt;
    stat.push_back(st);
}


//**************************** class DSagLearn ****************************

template<int stimd>
class DSagLearn {
public:
    using Agent = DSagent<stimd>;
    // using iarr = std::array<int, stimd>;
    using farr = Eigen::Matrix<flt, stimd, 1>;
    // using farr = std::array<flt, stimd>;
    using stat_type = LearnStat<stimd>;
    using vs_type = std::vector<stat_type>;
    using comp_stim = CompStim<stimd>;
    using rand_eng = std::mt19937;
    using rand_uni = std::uniform_real_distribution<flt>;
    using rand_int = std::uniform_int_distribution<int>;
    using rand_norm = std::normal_distribution<flt>;
    using rand_discr = std::discrete_distribution<int>;
    DSagLearn(int a_atyp,
              int a_repl,
              int a_T,
              flt a_sdR,
              flt a_sdx,
              flt a_avx1,
              flt a_avx2,
              flt a_w_true1,
              flt a_w_true2,
              flt a_w_true3,
              flt a_w_true4,
              flt a_w0,
              flt a_kapp0,
              flt a_s20,
              flt a_bet0,
              flt a_h0,
              flt a_s2R,
              flt a_omega,
              flt a_mu) :
        atyp{a_atyp},
        repl{a_repl},
        T{a_T},
        sdR{a_sdR},
        sdx{a_sdx},
        avx1{a_avx1},
        avx2{a_avx2},
        w_true1{a_w_true1},
        w_true2{a_w_true2},
        w_true3{a_w_true3},
        w_true4{a_w_true4},
        w0{a_w0},
        kapp0{a_kapp0},
        s20{a_s20},
        bet0{a_bet0},
        h0{a_h0},
        s2R{a_s2R},
        omega{a_omega},
        mu{a_mu},
        agent(w0, kapp0, s20, bet0, h0, s2R),
        cst1(sdR, 1),
        cst2(sdR, 2),
        cst3(sdR, 3),
        cst4(sdR, 4) { 
            stat.reserve(2*T);
            Setup_CompStim();
        }
    const vs_type& Get_stat() const { return stat; }
    void Setup_CompStim();
    void Learn(rand_eng& eng);
    void Add_stat(int tstep, const comp_stim& cs1, const comp_stim& cs2);
    int atyp;        // agent type
    int repl;        // replicate number (of agent type)
    int T;           // number of learning rounds
    flt sdR;         // sd of random reward variation
    flt sdx;         // sd of variation in x for first dimension
    flt avx1;        // mean of first dim of x for cs type 1
    flt avx2;        // mean of first dim of x for cs type 2
    flt w_true1;     // reward per x for 1st stimulus dimension
    flt w_true2;     // reward for presence for 2nd dimension
    flt w_true3;     // reward for presence for 3rd dimension
    flt w_true4;     // reward for presence for 3rd dimension
    flt w0;          // starting estimated feature value
    flt kapp0;       // starting feature learning rate
    flt s20;         // starting value for 'uncertainties in predictions'
    flt bet0;        // starting value for log(s2)
    flt h0;          // starting value for h quantity
    flt s2R;         // 'observation error' in the DS approach
    flt omega;       // parameter giving choice prob from values
    flt mu;          // meta-learning rate for DS approach
    Agent agent;     // learning agent
    comp_stim cst1;  // Type 1 compound stimulus
    comp_stim cst2;  // Type 2 compound stimulus
    comp_stim cst3;  // Type 3 compound stimulus
    comp_stim cst4;  // Type 4 compound stimulus
    vs_type stat;    // learning statistics
};

template<int stimd>
void DSagLearn<stimd>::Setup_CompStim()
{
    // set up compound stimulus types
    cst1.w_true[0] = w_true1;
    cst2.w_true[0] = w_true1;
    cst3.w_true[0] = w_true1;
    cst4.w_true[0] = w_true1;

    cst1.w_true[1] = w_true2;
    cst2.w_true[1] = w_true2;
    cst3.w_true[1] = w_true2;
    cst4.w_true[1] = w_true2;

    cst1.w_true[2] = w_true3;
    cst2.w_true[2] = w_true3;
    cst3.w_true[2] = w_true3;
    cst4.w_true[2] = w_true3;

    cst1.w_true[3] = w_true4;
    cst2.w_true[3] = w_true4;
    cst3.w_true[3] = w_true4;
    cst4.w_true[3] = w_true4;

    cst1.x[0] = avx1;
    cst2.x[0] = avx2;
    cst3.x[0] = avx2;
    cst4.x[0] = avx2;

    cst1.x[1] = 0;
    cst2.x[1] = 1;
    cst3.x[1] = 1;
    cst4.x[1] = 1;

    cst1.x[2] = 0;
    cst2.x[2] = 0;
    cst3.x[2] = 1;
    cst4.x[2] = 0;

    cst1.x[3] = 0;
    cst2.x[3] = 0;
    cst3.x[3] = 0;
    cst4.x[3] = 1;
}

template<int stimd>
void DSagLearn<stimd>::Learn(rand_eng& eng)
{
    // Simulate learning by Dayan-Sutton agent over T rounds. 
    rand_uni uni(0, 1);
    rand_norm nrmx(0, sdx); // random variation in first dimension
    rand_int int1_4(1, 4);  // uniform on {1, 2, 3, 4}

    // run through 2*T rounds of learning 
    for (int tstep = 1; tstep <= 2*T; ++tstep) {
        // set the two compound stimuli for this round (pedestrian!)
        comp_stim cs1 = cst1;
        comp_stim cs2 = cst2;
        // check if we are in first or second half of learning
        if (tstep <= T) {
            // random selection among first two types
            cs1 = (uni(eng) < 0.5) ? cst1 : cst2;
            cs2 = (uni(eng) < 0.5) ? cst1 : cst2;
        } else {
            // random selection among all 4 types
            int typ = int1_4(eng);
            if (typ == 1) cs1 = cst1;
            else if (typ == 2) cs1 = cst2;
            else if (typ == 3) cs1 = cst3;
            else cs1 = cst4;
            typ = int1_4(eng);
            if (typ == 1) cs2 = cst1;
            else if (typ == 2) cs2 = cst2;
            else if (typ == 3) cs2 = cst3;
            else cs2 = cst4;
        }
        // introduce random variation in x[0]
        cs1.x[0] += nrmx(eng);
        cs2.x[0] += nrmx(eng);
        // put the features states into the agent
        for (int i = 0; i < stimd; ++i) {
            agent.x1[i] = cs1.x[i];
            agent.x2[i] = cs2.x[i];
        }
        // compute estimated reward values
        agent.Q1 = agent.Qval(agent.x1);
        agent.Q2 = agent.Qval(agent.x2);
        // agent choice (based on estimated reward values)
        agent.choice = (uni(eng) < agent.PrCS1(omega)) ? 1 : 2;
        comp_stim& cs = (agent.choice == 1) ? cs1 : cs2;
        // get reward from selected compound stimulus
        agent.Rew = cs.Reward(eng);
        // store x and Q in agent
        agent.x = cs.x;
        agent.Q = (agent.choice == 1) ? agent.Q1 : agent.Q2;
        // prediction error
        agent.delt = agent.Rew - agent.Q;
        // update learning rates for the DS agent
        agent.DSupdate(mu);
        // update estimated feature values
        agent.Updatew();
        // store learning statistics for this round
        Add_stat(tstep, cs1, cs2);
    }
}

// Add data from current round to the learning statistics container 
template<int stimd>
void DSagLearn<stimd>::Add_stat(int tstep, const comp_stim& cs1, const comp_stim& cs2) 
{
    stat_type st;
    st.atyp = atyp;
    st.repl = repl;
    st.tstep = tstep;
    st.choice = agent.choice;
    st.cstyp1 = cs1.cstyp;
    st.cstyp2 = cs2.cstyp;
    st.x1 = agent.x1;
    st.x2 = agent.x2;
    st.x = agent.x;
    st.w_true = (st.choice == 1) ? cs1.w_true : cs2.w_true; 
    st.w = agent.w;
    st.kapp = agent.kapp;
    st.var1 = agent.s2;
    st.var2 = agent.bet;
    st.var3 = agent.h;
    st.Rew = agent.Rew;
    st.Q1 = agent.Q1;
    st.Q2 = agent.Q2;
    st.Q = agent.Q;
    st.Q1true = cs1.w_true.dot(cs1.x);
    st.Q2true = cs2.w_true.dot(cs2.x);
    st.delt = agent.delt;
    stat.push_back(st);
}

#endif // AGLEARN_HPP
