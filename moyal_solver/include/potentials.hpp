#pragma once
#include "moyal_setup.hpp"
#include<memory>
#include<cmath>

/**
 * @brief gaussian potential
 */
class GaussianPotential : public Potential {
    private:
        double amplitude_;
        double center_;
        double width_;

    public:
        GaussianPotential(double amplitude, double center, double width) : amplitude_(amplitude), center_(center), width_(width) {};

        double operator()(double x, double t = 0.0) const override {
            double dx = x - center_;
            return amplitude_ * std::exp(-dx*dx / (2.0 * width_*width_));
        }

        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<GaussianPotential>(amplitude_, center_, width_);
        }
};

/**
 * @brief harmonic oscillator potential
 * clone() in code, move-semantics
 */
class HarmonicPotential : public Potential {
    public:
    double mass_ = 1;
    double omega_ = 1;
    double center_ = 0.0;
        double operator()(double x, double t = 0.0) const override {
            double dx = x - center_;
            return 0.5 * mass_ * omega_ * omega_ * dx * dx;
        }

        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<HarmonicPotential>();
        }
};

/**
 * @brief Godbeer potential - 4th order polynomial
 */
class GodbeerAT : public Potential {
    public:
        double alpha = 1.963;
        double a4 = 0.0207;
        double a3 = -0.0053;
        double a2 = -0.0414;
        double a1 = 0.0158;
        double a0 = 0.0312;

        double operator()(double v, double t = 0.0) const override {
            double x = v/alpha;
            return a4*x * x * x * x + a3*x * x * x + a2*x * x + a1*x + a0;
        }
        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<GodbeerAT>();
        }
};

/**
 * @brief Godbeer potential - 4th order polynomial
 */
class SlocombeGC : public Potential {
    public:
        double V1=0.1617;
        double V2=0.082;
        double a1=0.305;
        double a2=0.755;
        double r1=-2.7;
        double r2=2.1;

        double const_term = 0.166 + 0.00019;

        double operator()(double q, double t = 0.0) const override {
            double exp1 = std::exp(-2 * a1 * (q - r1));
            double exp2 = std::exp(-a1 * (q - r1));
            double exp3 = std::exp(-2 * a2 * (r2 - q));
            double exp4 = std::exp(-a2 * (r2 - q));
            return V1 * (exp1 - 2 * exp2) + V2 * (exp3 - 2 * exp4) + const_term;
        }
        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<SlocombeGC>();
        }
};
/**
 * @brief free particle potential
 */
class FreeParticle : public Potential {
public:
    double operator()(double x, double t = 0.0) const override {
        return 0.0; 
    }
    
    std::unique_ptr<Potential> clone() const override {
        return std::make_unique<FreeParticle>();
    }
};

class ModelGC : public Potential {
    private:
    double a_can = 0.006457467585167605;
    double b_can = 0.03131130708309966;
    double c_can = 0.03979932858576989;

    double a_bar = -0.006425438605566449;
    double b_bar = 0.0038910266591298437;
    double c_bar = 0.025223914926830745;

    double a_tau = 0.013406834311699567;
    double b_tau = -0.044575846347551726;
    double c_tau = 0.05473263795143025;

    public:
        ModelGC() {};

        double operator()(double x, double t = 0.0) const override {
            if (x < -1.03){
                return a_can * x * x + b_can * x + c_can;
            }
            else if (x < 1.15){
                return a_bar * x * x + b_bar * x + c_bar;
            }
            else {
                return a_tau * x * x + b_tau * x + c_tau;
            } 
        }

        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<ModelGC>();
        }
};

class ModelAT : public Potential {
    private:
    
    double a_can = 0.01548757014342916;
    double b_can = 0.0544195348802887;
    double c_can = 0.04817678375503837;

    double a_bar = -0.010377758242695038;
    double b_bar = 0.007613726939775605;
    double c_bar = 0.0310333936186076;

    double a_tau = 0.0125386460155263;
    double b_tau = -0.04565578211748337;
    double c_tau = 0.0619375921034291;

    public:
        ModelAT() {};

        double operator()(double x, double t = 0.0) const override {
            if (x < -0.51){
                return a_can * x * x + b_can * x + c_can;
            }
            else if (x < 1.21){
                return a_bar * x * x + b_bar * x + c_bar;
            }
            else {
                return a_tau * x * x + b_tau * x + c_tau;
            } 
        }

        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<ModelAT>();
        }
};
