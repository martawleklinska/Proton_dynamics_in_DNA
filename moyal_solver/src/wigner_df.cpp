#include "wigner_df.hpp"
#include <cmath>
#include <fstream>
#include<iostream>
#include <stdexcept>

WDF::WDF(const PhaseSpace& phase_space) 
    : phase_space_(phase_space) {
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    wigner_.resize(nx, std::vector<Complex>(np, Complex(0, 0)));
}

void WDF::initializeGaussian(double x0, double p0, double sigma_x, double sigma_p) {
    const auto& X = phase_space_.X();
    const auto& P = phase_space_.P();
    const double hbar = 1.0; // From config
    
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            double x = X[i][j] - x0;
            double p = P[i][j] - p0;
            
            double exponent = -x*x / (2.0 * sigma_x*sigma_x) 
                             - 2.0 * sigma_x*sigma_x * p*p / (hbar*hbar);
            
            // Note: ignoring sigma_p so that to match \delta x\delta p \sim \hbar
            (void)sigma_p; // Suppress warning
            
            wigner_[i][j] = std::exp(exponent) / (M_PI * hbar);
        }
    }
    normalize();
}

void WDF::initializeCoherentState(double x0, double p0, double sigma) {
    // Coherent state is a special case of Gaussian with sigma_x = sigma, sigma_p = ħ/(2σ)
    double sigma_p = 1.0 / (2.0 * sigma); // Assuming ħ = 1
    initializeGaussian(x0, p0, sigma, sigma_p);
}

void WDF::initializeFromFunction(const std::function<Complex(double, double)>& func) {
    const auto& X = phase_space_.X();
    const auto& P = phase_space_.P();
    
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            wigner_[i][j] = func(X[i][j], P[i][j]);
        }
    }
}

double WDF::norm() const {
    double sum = 0.0;
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    double dx = phase_space_.dx();
    double dp = phase_space_.dp();
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            sum += std::abs(wigner_[i][j]) * dx * dp;
        }
    }
    return sum;
}

double WDF::totalProbability() const {
    double sum = 0.0;
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    double dx = phase_space_.dx();
    double dp = phase_space_.dp();
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            sum += wigner_[i][j].real() * dx * dp;  
        }
    }
    return sum;
}
std::pair<double, double> WDF::expectedPosition() const {
    const auto& X = phase_space_.X();
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    double dx = phase_space_.dx();
    double dp = phase_space_.dp();
    
    double mean_x = 0.0;
    double total_prob = 0.0;
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            double w_real = wigner_[i][j].real();
            if (w_real > 0) {
                mean_x += X[i][j] * w_real * dx * dp;
                total_prob += w_real * dx * dp;
            }
        }
    }
    
    if (total_prob > 0) {
        mean_x /= total_prob;
    }
    
    return {mean_x, total_prob};
}

std::pair<double, double> WDF::expectedMomentum() const {
    const auto& P = phase_space_.P();
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    double dx = phase_space_.dx();
    double dp = phase_space_.dp();
    
    double mean_p = 0.0;
    double total_prob = 0.0;
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            double w_real = wigner_[i][j].real();
            if (w_real > 0) {
                mean_p += P[i][j] * w_real * dx * dp;
                total_prob += w_real * dx * dp;
            }
        }
    }
    
    if (total_prob > 0) {
        mean_p /= total_prob;
    }
    
    return {mean_p, total_prob};
}

void WDF::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    const auto& X = phase_space_.X();
    const auto& P = phase_space_.P();
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            file << X[i][j] << " " << P[i][j] << " " 
                 << wigner_[i][j].real() << " " << wigner_[i][j].imag() << "\n";
        }
        file << "\n";
    }
}

RealMatrix WDF::real() const {
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    RealMatrix result(nx, std::vector<double>(np));
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            result[i][j] = wigner_[i][j].real();
        }
    }
    return result;
}

RealMatrix WDF::magnitude() const {
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    RealMatrix result(nx, std::vector<double>(np));
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            result[i][j] = std::abs(wigner_[i][j]);
        }
    }
    return result;
}

void WDF::normalize() {
    double total_prob = totalProbability();  // <-- użyj totalProbability(), nie norm()!
    
    if (std::abs(total_prob) > 1e-15) {  // sprawdź czy nie zero
        int nx = phase_space_.gridX();
        int np = phase_space_.gridP();
        
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < np; ++j) {
                wigner_[i][j] /= total_prob;
            }
        }
        
        std::cout << "Normalized: initial probability was " << total_prob << "\n";
    } else {
        std::cerr << "ERROR: Cannot normalize - total probability is zero!\n";
    }
}
