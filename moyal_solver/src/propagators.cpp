#include"propagators.hpp"
#include<cmath>

KineticPropagator::KineticPropagator(const PhaseSpace& phase_space, double mass)
    : Propagator(phase_space), mass_(mass) {
}

void KineticPropagator::computePropagator(double dt){
    const auto& P = phase_space_.P();  

    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();

    for (int i = 0; i < nx; i++){
        for (int j = 0; j < np; j++){
            double p = P[i][j];
            double kinetic_energy = p * p / (2.0 * mass_);
            propagator_matrix_[i][j] = std::exp(Complex(0, -dt * kinetic_energy));
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