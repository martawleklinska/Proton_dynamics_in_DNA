#include "moyal_solver.hpp"
#include "potentials.hpp"
#include <iostream>
#include <filesystem>

int main() {
    try {
        // Create output directory
        std::filesystem::create_directories("./output");
        
        // Configuration 
        MoyalConfig config;
        
        std::cout << "=== Moyal Equation Solver ===\n";
        std::cout << "Grid: " << config.gridX << " x " << config.gridP << "\n";
        std::cout << "Time steps: " << config.timeSteps << ", dt = " << config.dt << "\n";
        std::cout << "Initial conditions: x0 = " << config.x_init 
                  << ", p0 = " << config.p_init << "\n\n";
        
        // Create potential - we're choosing Duffing
        auto potential = std::make_unique<GodbeerAT>();
        
        // Create solver
        MoyalSolver solver(config, std::move(potential));
        
        // Set initial Gaussian state
        solver.setInitialGaussian(config.x_init, config.p_init, 
                                 config.sigma_x, config.sigma_p);
        
        std::cout << "Initial state set. Starting evolution...\n";
        
        // Run evolution
        solver.evolve();
        
        std::cout << "\nEvolution completed successfully!\n";
        std::cout << "Output files saved to ./output/\n";
        
        // Show final statistics
        double mean_x, mean_p, sigma_x, sigma_p;
        solver.computeExpectationValues(mean_x, mean_p, sigma_x, sigma_p);
        
        std::cout << "\nFinal state statistics:\n";
        std::cout << "Mean position: " << mean_x << "\n";
        std::cout << "Mean momentum: " << mean_p << "\n";
        std::cout << "Position uncertainty: " << sigma_x << "\n";
        std::cout << "Momentum uncertainty: " << sigma_p << "\n";
        std::cout << "Total probability: " << solver.getWigner().totalProbability() << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
