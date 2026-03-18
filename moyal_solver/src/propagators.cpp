#include"propagators.hpp"
#include<cmath>

KineticPropagator::KineticPropagator(const PhaseSpace& phase_space, double mass)
    : Propagator(phase_space), mass_(mass) {
}

void KineticPropagator::computePropagator(double dt) {
    // Po fft_x, dane są w przestrzeni (lambda, p)
    // lambda = KX, p = P
    // Propagator: exp(-i * p * lambda / mass * dt)
    
    const auto& KX = phase_space_.KX();  // lambda - sprzężone do x
    const auto& P  = phase_space_.P();   // rzeczywisty pęd
    
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < np; j++) {
            double lambda = KX[i][j];
            double p      = P[i][j];
            // dW/dt = -(p/m) * d/dx W  =>  w przestrzeni Fouriera: -i*(p/m)*lambda
            double phase = -dt * p * lambda / mass_;
            propagator_matrix_[i][j] = std::exp(Complex(0.0, phase));
        }
    }
}

void KineticPropagator::apply(ComplexMatrix& wigner, double dt){
    computePropagator(dt);

    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();

    for (int i = 0; i < nx; i++){
        for (int j = 0; j < np; j++){
            wigner[i][j] *= propagator_matrix_[i][j];
        }
    }
}

PotentialPropagator::PotentialPropagator(const PhaseSpace& phase_space,
                    std::unique_ptr<Potential> potential,
                    double hbar)
    : Propagator(phase_space), potential_(std::move(potential)), hbar_(hbar){

    }

void PotentialPropagator::computePropagator(double /* dt*/){

}
void PotentialPropagator::apply(ComplexMatrix& wigner, double dt){
    computeMoyalBracket(wigner, dt);
}
void PotentialPropagator::computeMoyalBracket(ComplexMatrix& wigner, double dt, double t){
    const auto& X = phase_space_.X();
    const auto& KP = phase_space_.KP();

    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();

    for (int i = 0; i < nx; i++){
        for (int j = 0; j < np; j++){
            double x = X[i][j];
            double theta = 0.5 * hbar_ * KP[i][j];

            double V_plus = (*potential_)(x + theta, t);
            double V_minus = (*potential_)(x - theta, t);
            double potential_diff = V_plus - V_minus;

            Complex phase_factor = std::exp(Complex(0, -dt * potential_diff / hbar_));
            propagator_matrix_[i][j] = phase_factor;
            wigner[i][j] *=phase_factor;
        }
    }
}