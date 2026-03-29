#include "moyal_solver.hpp"
#include "potentials.hpp"
#include <iostream>
#include <filesystem>

int main() {

        std::filesystem::create_directories("./output");
        
        MoyalConfig config;
        
        std::cout << "=== Moyal Equation Solver ===\n";
        std::cout << "Grid: " << config.gridX << " x " << config.gridP << "\n";
        std::cout << "Time steps: " << config.timeSteps << ", dt = " << config.dt << "\n";
        std::cout << "Initial conditions: x0 = " << config.x_init 
                  << ", p0 = " << config.p_init << "\n\n";
        
        auto potential = std::make_unique<HarmonicPotential>();
        
        MoyalSolver solver(config, std::move(potential));
        
        solver.setInitialGaussian(config.x_init, config.p_init, 
                                 config.sigma_x, config.sigma_p);
        
        std::cout << "Initial state set. Starting evolution...\n";
        
        solver.evolve();
        

    
    return 0;
}
