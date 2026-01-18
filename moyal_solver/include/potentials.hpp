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
    double center_ = 0.5;
        double operator()(double x, double t = 0.0) const override {
            double dx = x - center_;
            return 0.5 * mass_ * omega_ * omega_ * dx * dx;
        }

        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<HarmonicPotential>();
        }
};

/**
 * @brief DK potential exp
 */
class OGPotential : public Potential {
    public:
        double operator()(double x, double t = 0.0) const override {
            return 0.008 * std::exp(-std::pow(x + 200.0, 2) / 250.);
        }
        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<OGPotential>();
        }
};
/**
 * @brief Duffing potential
 */
class Duffing : public Potential {
public:
    double alpha = -5e-1;
    double beta  =  1e-3;
    double gamma =  0.5;
    double omega =  1.0e-1;


    double operator()(double x, double t = 0.0) const override {
        return 0.5*alpha*x*x + 0.25*beta*x*x*x*x - x*gamma * std::cos(omega*t);  
    }
    
    std::unique_ptr<Potential> clone() const override {
        return std::make_unique<Duffing>();
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
            return std::make_unique<GodbeerAT>();
        }
};