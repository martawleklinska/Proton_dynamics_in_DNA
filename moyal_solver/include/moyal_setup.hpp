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

    int gridX = 1024;
    int gridP = 1024;
    double ampX = 30.0;   
    double ampP = 100.0;   
    
    double hbar = 1.0;
    double mass = 1.0;
    
    double dt = 0.01;      
    int timeSteps = 900;
    
    double x_init  = -1.0;
    double p_init  = 1.0;
    double sigma_x = 0.5;
    double sigma_p = 1.0;
    
    int outputEvery = 10;
    std::string outputDir = "./output/";
};