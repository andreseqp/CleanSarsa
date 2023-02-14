#ifndef AGENT_HPP
#define AGENT_HPP

// #include "NNetw.hpp"
#include <string>
#include <array>
#include <ostream>
#include <istream>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>

using flt = double;

// Concerning learning rates for the agents, the code below uses a notation
// where kapp is a learning rate that includes the multiplicative contribution
// from the stimulus value (the cue) x, whereas for the Rescorla-Wagner the
// notation alph is used for the learning rate without the contribution from x
// (this seems to be a 'standard convention' in the literature)


//************************* struct RWagent ******************************

// Agent with Rescorla-Wagner learning

template<int stimd>
struct RWagent {
public:
    using farr = Eigen::Matrix<flt, stimd, 1>;
    RWagent(flt w0 = 0,
            flt alph0 = 1) //:
    {
        w.fill(w0);
        alph.fill(alph0);
    }
    // return the estimated reward value of feature vector x
    flt Qval(const farr& xv) const 
    { 
        return w.dot(xv);
    }
    // return the probability of selecting compound stimulus 1 (sigmoid)
    flt PrCS1(flt omega) const { return 1/(1 + std::exp(-omega*(Q1 - Q2))); }
    // Learning update of w, using current learning rates
    void Updatew();
    // public data members
    farr w;      // estimates of reward values for features
    farr alph;   // learning rate
    farr x1;     // features of compound stimulus 1
    farr x2;     // features of compound stimulus 2
    farr x;      // features of selected compound stimulus
    flt Rew;     // perceived reward
    flt Q1;      // estimated value compound stimulus 1
    flt Q2;      // estimated value compound stimulus 2
    flt Q;       // estimated value of selected compound stimulus
    flt delt;    // 'prediction error' (TD error)
    int choice;  // choice of compound stimulus (1 or 2)
};

template<int stimd>
void RWagent<stimd>::Updatew()
{
    for (int i = 0; i < stimd; ++i) {
        w[i] += alph[i]*x[i]*delt;
    }
}


//************************* struct PMagent ******************************

// Agent with Pearce-Macintosh learning

template<int stimd>
struct PMagent {
public:
    using farr = Eigen::Matrix<flt, stimd, 1>;
    // constructor
    PMagent(flt w0 = 0,
            flt kapp0 = 1,
            flt alph0 = 1,
            flt sigm0 = 1,
            flt beta0 = 1) :
        beta{beta0}
    {
        wprev.fill(w0);
        w.fill(w0);
        kapp.fill(kapp0);
        alph.fill(alph0);
        sigm.fill(sigm0);
    }
    // return the estimated reward value of feature vector x
    flt Qval(const farr& xv) const 
    { 
        return w.dot(xv);
    }
    // return the probability of selecting compound stimulus 1 (sigmoid)
    flt PrCS1(flt omega) const { return 1/(1 + std::exp(-omega*(Q1 - Q2))); }
    // Learning update of w, using current learning rates
    void Updatew();
    // Pearce-Mackintosh 'unified' meta learning update
    void PMupdate(flt gam_a, flt gam_s);
    // public data members
    farr wprev;  // previous estimates of reward values for features
    farr w;      // estimates of reward values for features
    farr kapp;   // learning rate
    farr alph;   // learning rate
    farr sigm;   // learning rate
    farr x1;     // features (0/1) of compound stimulus 1
    farr x2;     // features (0/1) of compound stimulus 2
    farr x;      // features (0/1) of selected compound stimulus
    flt Rew;     // perceived reward
    flt Q1;      // estimated value compound stimulus 1
    flt Q2;      // estimated value compound stimulus 2
    flt Q;       // estimated value of selected compound stimulus
    flt delt;    // 'prediction error' (TD error)
    flt beta;    // learning rate
    int choice;  // choice of compound stimulus (1 or 2)
};

template<int stimd>
void PMagent<stimd>::Updatew()
{
    wprev = w;
    for (int i = 0; i < stimd; ++i) {
        w[i] += alph[i]*sigm[i]*beta*x[i]*delt;
    }
}

template<int stimd>
void PMagent<stimd>::PMupdate(flt gam_a, flt gam_s)
{
    for (int i = 0; i < stimd; ++i) {
        flt adiff1 = std::abs(Rew - wprev[i]*x[i]);
        flt adiff2 = std::abs(Rew - (Q - wprev[i]*x[i]));
        // flt adiff1 = std::abs(Rew - w[i]*x[i]);
        // flt adiff2 = std::abs(Rew - (Q - w[i]*x[i]));

        alph[i] += gam_a*(adiff2 - adiff1);
        sigm[i] = (1 - gam_s)*sigm[i] + gam_s*adiff1;

        // impose lower and upper limits on alph and sigm, following Pearce
        // and Mackintosh (2010)
        if (alph[i] > 1) {
            alph[i] = 1;
        } else if (alph[i] < 0.05) {
            alph[i] = 0.05;
        }
        if (sigm[i] > 1) {
            sigm[i] = 1;
        } else if (sigm[i] < 0.5) {
            sigm[i] = 0.5;
        }
        kapp[i] = alph[i]*sigm[i]*beta*x[i];
    }
}


//************************* struct DSagent ******************************

// Agent with Dayan-Sutton (Kalman filter inspired) learning

template<int stimd>
struct DSagent {
public:
    using farr = Eigen::Matrix<flt, stimd, 1>;
    using fmat = Eigen::Matrix<flt, stimd, stimd>;
    // constructor
    DSagent(flt w0 = 0,
            flt kapp0 = 1,
            flt a_s20 = 1,
            flt bet0 = 0,
            flt h0 = 0,
            flt a_s2R = 1,
            flt a_tau2 = 0) :
        s20{a_s20},
        s2R{a_s2R},
        tau2{a_tau2}
    {
        w.fill(w0);
        kapp.fill(kapp0);
        s2.fill(s20);
        bet.fill(bet0);
        h.fill(h0);
        Sig.fill(0);
        for (int i = 0; i < stimd; ++i) {
            Sig(i, i) = s2[i];
        }
    }
    // return the estimated reward value of feature vector x
    flt Qval(const farr& xv) const 
    { 
        return w.dot(xv);
    }
    // return the probability of selecting compound stimulus 1 (sigmoid)
    flt PrCS1(flt omega) const { return 1/(1 + std::exp(-omega*(Q1 - Q2))); }
    // Learning update of w, using current learning rates
    void Updatew();
    // ReLU function
    flt ReLU(flt X) const { return (X > 0) ? X : 0; }
    // Dayan-Sutton (Kalman filter inspired) meta-learning update
    void DSupdate(flt mu);
    // public data members
    farr w;      // estimates of reward values for features
    farr kapp;   // learning rate
    farr s2;     // 'uncertainties in predictions'
    farr bet;    // log(s2)
    farr h;      // quantity introduced by Sutton (1992)
    farr x1;     // feature states of compound stimulus 1
    farr x2;     // feature states of compound stimulus 2
    farr x;      // feature states of selected compound stimulus
    fmat Sig;    // Variance covariance matrix (Kalman filter)
    flt Rew;     // perceived reward
    flt Q1;      // estimated value compound stimulus 1
    flt Q2;      // estimated value compound stimulus 2
    flt Q;       // estimated value of selected compound stimulus
    flt delt;    // 'prediction error' (TD error)
    flt s20;     // starting values for s2
    flt s2R;     // 'observation error' in the DS algorithm
    flt tau2;    // variance of increments to w
    int choice;  // choice of compound stimulus (1 or 2)
};

template<int stimd>
void DSagent<stimd>::Updatew()
{
    for (int i = 0; i < stimd; ++i) {
        w[i] += kapp[i]*delt;
    }
}

template<int stimd>
void DSagent<stimd>::DSupdate(flt mu)
{
    // this is the K1 algorithm of Sutton 1992b
    flt sums2x2 = 0;
    for (int i = 0; i < stimd; ++i) {
        bet[i] += mu*delt*h[i]*x[i];  // equation (13) in Sutton 1992b
        s2[i] = std::exp(bet[i]);     // equation (11) in Sutton 1992b
        sums2x2 += s2[i]*x[i]*x[i];
    }
    for (int i = 0; i < stimd; ++i) {
        kapp[i] = s2[i]*x[i]/(sums2x2 + s2R); // equation (10) in Sutton 1992b
        h[i] = (h[i] + kapp[i]*delt)*ReLU(1 - kapp[i]*x[i]); // equation (15)
    }

    // // this is the 'original' Kalman algorithm (from Sutton 1992b or Gershman 2015);
    // // first add effect of random increment to Sig
    // fmat Incr = fmat::Zero();
    // for (int i = 0; i < 5; ++i) {
    //     Incr(i, i) = tau2;
    // }
    // Sig = Sig + Incr; // Add effect of random increment to Sig
    // farr Sigx = Sig*x; // numerator in expression for kapp
    // flt den = x.transpose()*Sigx + s2R; // denominator
    // kapp = Sigx/den; // this is the Kalman gain (learning rates)
    // // update variance-covariance matrix
    // Sig = (fmat::Identity() - kapp*x.transpose())*Sig;
    // // finally transfer diagonal elements to s2
    // for (int i = 0; i < stimd; ++i) {
    //     s2[i] = Sig(i, i);
    // }

    // // this is the IDBD algorithm of Sutton 1992a, as described in Sutton
    // // 1992b (there is a difference in notation, involving whether the cue x
    // // is part of the learning rate); note that starting value of bet needs to
    // // be right to prevent divergence of the algorithm (bet0 = log(1/stimd)
    // // according to Sutton 1992b)
    // for (int i = 0; i < stimd; ++i) {
    //     bet[i] += mu*delt*h[i]*x[i];
    //     kapp[i] = std::exp(bet[i])*x[i];  // eq. (17) in Sutton 1992b, cf (4) in Sutton 1992a
    // }
    // for (int i = 0; i < stimd; ++i) {
    //     // from eq. (20) in Sutton 1992b, cf (5) or fig 2 in Sutton 1992a
    //     h[i] = h[i]*ReLU(1 - kapp[i]*x[i]) + kapp[i]*delt;
    // }

    // // this is the Yu and Dayan algorithm (from Yu and Dayan 1993 and Dayan and Yu 1993);
    // // Assume there is no random increment to Sig (stable situation between change points)
    // // fmat Incr = fmat::Zero();
    // // for (int i = 0; i < ndim; ++i) {
    // //     Incr(i, i) = tau2; // tau2 is zero by default
    // // }
    // // Sig = Sig + Incr; // Add effect of random increment to Sig
    // farr Sigx = Sig*x; // numerator in expression for kapp
    // flt xSigx = x.transpose()*Sigx; 
    // flt den = xSigx + s2R; // denominator
    // kapp = Sigx/den; // this is the Kalman gain (learning rates)
    // // update variance-covariance matrix
    // Sig = (fmat::Identity() - kapp*x.transpose())*Sig;
    // // finally transfer diagonal elements to s2
    // for (int i = 0; i < stimd; ++i) {
    //     s2[i] = Sig(i, i);
    // }
    // // check is there is a change point
    // if (delt*delt/xSigx > 3) {
    //     // reset Sig to 'starting value'
    //     Sig.fill(0);
    //     for (int i = 0; i < stimd; ++i) {
    //         Sig(i, i) = s20;
    //     }
    // }
}



// //************************* struct NNagent ******************************

// template<int stimd, int inpd, int hiddend, int outpd>
// struct NNagent {
// public:
//     using Netw = NN0<inpd, hiddend, outpd>;
//     using flt = typename Netw::flt;
//     using VecIn = typename Netw::VecIn;
//     using VecHid = typename Netw::VecHid;
//     using VecOut = typename Netw::VecOut;
//     using Mat1 = typename Netw::Mat1;
//     using Mat2 = typename Netw::Mat2;
//     using iVec = Eigen::Matrix<int, stimd, 1>;
//     using farr = Eigen::Matrix<flt, stimd, 1>;
//     // constructor
//     NNagent(const VecIn& a_x, const Mat1& a_m1, const Mat2& a_m2) :
//         ntw(a_x, a_m1, a_m2) {}
//     // public data members
//     Netw ntw;    // network for decision making
//     int choice;  // choice of compound stimulus (1 or 2)
// };


#endif // AGENT_HPP
