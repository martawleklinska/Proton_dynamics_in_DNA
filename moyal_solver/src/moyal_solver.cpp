#include "moyal_solver.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>

MoyalSolver::MoyalSolver(const MoyalConfig& config, std::unique_ptr<Potential> potential)
    : config_(config), phase_space_(config), wigner_(phase_space_),
      current_step_(0), current_time_(0.0) {
    
    kinetic_prop_ = std::make_unique<KineticPropagator>(phase_space_, config_.mass);
    potential_prop_ = std::make_unique<PotentialPropagator>(phase_space_, 
                                                           std::move(potential), 
                                                           config_.hbar);
}
// move semantics
void MoyalSolver::setInitialCondition(const std::function<Complex(double, double)>& init_func) {
    wigner_.initializeFromFunction(init_func);
}

void MoyalSolver::setInitialGaussian(double x0, double p0, double sigma_x, double sigma_p) {
    wigner_.initializeGaussian(x0, p0, sigma_x, sigma_p);
}

void MoyalSolver::evolve(int steps) {
    int total_steps = (steps > 0) ? steps : config_.timeSteps;
    
    std::cout << "Starting evolution for " << total_steps << " steps...\n";
    
    // Save initial state
    if (current_step_ % config_.outputEvery == 0) {
        outputData();
    }
    
    for (int step = 0; step < total_steps; ++step) {
        evolveOneStep();
        
        if ((current_step_) % config_.outputEvery == 0) {
            outputData();
        }
        
        if (current_step_ % 100 == 0) {
            std::cout << "Step " << current_step_ << "/" << total_steps 
                     << ", time = " << current_time_ << std::endl;
        }
    }
    
    std::cout << "Evolution completed.\n";
}

void MoyalSolver::evolveOneStep() {
    strangSplittingStep();
    current_step_++;
    current_time_ += config_.dt;
}

void MoyalSolver::strangSplittingStep() {
    // Proper Moyal equation Strang splitting:
    // exp(-iΔt/(2ℏ)T) exp(-iΔt/ℏ V) exp(-iΔt/(2ℏ)T)
    // Where T acts in momentum space, V acts in position space
    
    auto& wigner_data = wigner_.data();
    
    // 1. Half kinetic step: T(dt/2) in momentum space
    phase_space_.fft_x(wigner_data, true);   // FFT x to lambda (momentum space)
    kinetic_prop_->apply(wigner_data, config_.dt / 2.0);  // Apply T(dt/2)
    phase_space_.fft_x(wigner_data, false);  // IFFT lambda to x
    
    // 2. Full potential step: V(dt) in position space  
    phase_space_.fft_p(wigner_data, true);   // FFT p → θ (position space)
    potential_prop_->apply(wigner_data, config_.dt);      // Apply V(dt)
    phase_space_.fft_p(wigner_data, false);  // IFFT θ → p
    
    // 3. Half kinetic step: T(dt/2) in momentum space
    phase_space_.fft_x(wigner_data, true);   // FFT x → λ (momentum space)
    kinetic_prop_->apply(wigner_data, config_.dt / 2.0);  // Apply T(dt/2)
    phase_space_.fft_x(wigner_data, false);  // IFFT λ → x
    
}

void MoyalSolver::computeExpectationValues(double& mean_x, double& mean_p, 
                                         double& sigma_x, double& sigma_p) const {
    auto pos_result = wigner_.expectedPosition();
    auto mom_result = wigner_.expectedMomentum();
    
    mean_x = pos_result.first;
    mean_p = mom_result.first;
    
    // Calculate uncertainties 
    const auto& X = phase_space_.X();
    const auto& P = phase_space_.P();
    const auto& wigner_data = wigner_.data();
    
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    double dx = phase_space_.dx();
    double dp = phase_space_.dp();
    
    double var_x = 0.0, var_p = 0.0, total_prob = 0.0;
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            double w_real = wigner_data[i][j].real();
            if (w_real > 0) {
                double dx_dev = X[i][j] - mean_x;
                double dp_dev = P[i][j] - mean_p;
                
                var_x += dx_dev * dx_dev * w_real * dx * dp;
                var_p += dp_dev * dp_dev * w_real * dx * dp;
                total_prob += w_real * dx * dp;
            }
        }
    }
    
    if (total_prob > 0) {
        sigma_x = std::sqrt(var_x / total_prob);
        sigma_p = std::sqrt(var_p / total_prob);
    } else {
        sigma_x = sigma_p = 0.0;
    }
}

void MoyalSolver::saveState(const std::string& filename) const {
    wigner_.saveToFile(filename);
}

void MoyalSolver::outputData() {
    // Save current Wigner distribution
    std::ostringstream filename;
    filename << config_.outputDir << "wigner_" << std::setw(11) << std::setfill('0') 
             << current_step_ << ".dat";
    
    saveState(filename.str());
    
    // Calculate and save expectation values and nonclassicality parameter
    double mean_x, mean_p, sigma_x, sigma_p;
    computeExpectationValues(mean_x, mean_p, sigma_x, sigma_p);
    double delta = calculateNonclassicalityParameter();
    
    std::ofstream stats_file(config_.outputDir + "stats.dat", std::ios::app);
    if (current_step_ == 0) {
        std::cout << "data saved:# Step Time Mean_X Mean_P Sigma_X Sigma_P Norm Delta\n";
    }
    
    stats_file << current_step_ << " " << current_time_ << " "
               << mean_x << " " << mean_p << " " 
               << sigma_x << " " << sigma_p << " "
               << wigner_.totalProbability() << " " 
               << delta << "\n";

}

double MoyalSolver::calculateNonclassicalityParameter() const {
    const auto& wigner_data = wigner_.data();
    const int nx = phase_space_.gridX();
    const int np = phase_space_.gridP();
    const double dx = phase_space_.dx();
    const double dp = phase_space_.dp();
    
    double integral = 0.0;
    
    // Numerical integration using trapezoidal rule
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            // Calculate |W(x_i, p_j)|
            double wigner_abs = std::abs(wigner_data[i][j]);
            
            // Integration weights for trapezoidal rule
            double weight_x = (i == 0 || i == nx - 1) ? 0.5 : 1.0;
            double weight_p = (j == 0 || j == np - 1) ? 0.5 : 1.0;
            
            integral += weight_x * weight_p * wigner_abs * dx * dp;
        }
    }
    
    // \delta(t) = \iint |W(x,p,t)| dx dp - 1
    return integral - 1.0;
}
