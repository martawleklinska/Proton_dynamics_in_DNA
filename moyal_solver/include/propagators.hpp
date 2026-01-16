#pragma once
#include "moyal_setup.hpp"
#include "phase_space.hpp"

/**
 * @brief Base class for evolution propagators in Strang split-step method
 * Abstract class 
 */
class Propagator {
    protected:
        const PhaseSpace& phase_space_; ///< Reference to phase space grid
        ComplexMatrix propagator_matrix_; ///< Pre-computed propagation matrix
    
    public:
        /**
         * @brief Constructor
         * @param phase_space Reference to the phase space configuration
         */
        explicit Propagator(const PhaseSpace& phase_space) : phase_space_(phase_space) {
            int nx = phase_space_.gridX();
            int np = phase_space_.gridP();
            propagator_matrix_.resize(nx, std::vector<Complex>(np));
        }

        virtual ~Propagator() = default;

        /**
         * @brief Apply propagator to Wigner distribution function
         * @param wigner The Wigner function to evolve (modified in-place)
         * @param dt Time step for evolution
         */
        virtual void apply(ComplexMatrix& wigner, double dt) = 0;

        /**
         * @brief Access to pre-computed propagation matrix
         * @return Const reference to the propagator matrix
         */
        const ComplexMatrix& matrix() const {return propagator_matrix_;}
    
    protected:
        /**
         * @brief Compute the propagation matrix for given time step
         * @param dt Time step for which to compute the propagator
         */
        virtual void computePropagator(double dt) = 0;

};

/**
 * @brief Kinetic energy propagator for free particle evolution
 */
class KineticPropagator : public Propagator {
    private:
        double mass_; ///< Particle mass

    public:
        /**
         * @brief Constructor
         * @param phase_space Reference to phase space configuration
         * @param mass Particle mass for kinetic energy calculation
         */
        KineticPropagator(const PhaseSpace& phase_space, double mass);

        void apply(ComplexMatrix& wigner, double dt) override;

    protected:
        /**
         * @brief Compute kinetic propagator matrix exp(-iT*dt/ℏ)
         * @param dt Time step for propagator computation
         */
        void computePropagator(double dt) override;
};

/**
 * @brief Potential energy propagator using Moyal bracket formalism
 */
class PotentialPropagator : public Propagator {
    private:
        std::unique_ptr<Potential> potential_; ///< Potential energy function
        double hbar_; ///< Reduced Planck constant

    public:
        /**
         * @brief Constructor
         * @param phase_space Reference to phase space configuration
         * @param potential Unique pointer to potential energy function
         * @param hbar Reduced Planck constant 
         */
        PotentialPropagator(const PhaseSpace& phase_space,
                            std::unique_ptr<Potential> potential,
                            double hbar = 1.0);
        
        /**
         * @brief Apply potential evolution step to Wigner function
         * @param wigner Wigner function to evolve (modified in-place)
         * @param dt Time step for evolution
         */
        void apply(ComplexMatrix& wigner, double dt) override;

    protected:
        /**
         * @brief Compute potential propagator using Moyal bracket
         * @param dt Time step for propagator computation
         */
        void computePropagator(double dt) override;

    private:
        /**
         * @brief Compute Moyal bracket {H,W} for potential evolution
         * @param wigner Wigner function (modified in-place)
         * @param dt Time step for evolution
         * @param t Current time for time-dependent potentials (default: 0.0)
         */
        void computeMoyalBracket(ComplexMatrix& wigner, double dt, double t = 0.0);
};