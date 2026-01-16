#pragma once
#include<vector>
#include<complex>
#include<functional>
#include<memory>

using Complex = std::complex<double>;
using ComplexMatrix = std::vector<std::vector<Complex>>;
using RealMatrix = std::vector<std::vector<double>>;

/**
 * @brief Abstract Potential functions class
 * destructor must be virtual
 */
class Potential {
    public:
        virtual ~Potential() = default;
        virtual double operator()(double x, double t = 0.0) const = 0;
        virtual std::unique_ptr<Potential> clone() const = 0;
};

/**
 * @brief config params for Moyal eqn solver
 */

struct MoyalConfig {
    // grid params
    int gridX = 2048;      
    int gridP = 2048;
    double ampX = 50.0;   
    double ampP = 50.0;

    double hbar = 1.0;
    double mass = 1.0;

    // time params
    double dt = 0.1;
    int timeSteps = 200;

    // init conditions
    double x_init  = -2.;
    double p_init  = 1.15;

    double sigma_x = 1.0; 
    double sigma_p = 1/2;    ///< sigma_p = hbar/(2*sigma_x) = 1/(2*1.0)

    // output
    int outputEvery = 1;
    std::string outputDir = "./output/";
};