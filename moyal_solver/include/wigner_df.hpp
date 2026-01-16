#pragma once
#include"moyal_setup.hpp"
#include"phase_space.hpp"

/**
 * @brief WDF class
 */
class WDF{
    private:
        const PhaseSpace& phase_space_; ///< Reference to phase space grid
        ComplexMatrix wigner_;          ///< Pre-computed propagation matrix

    public:
        /**
         * @brief Constructor
         * @param phase_space Reference to the phase space configuration
         */
        explicit WDF(const PhaseSpace& phase_space);

        ComplexMatrix& data() {return wigner_;}
        const ComplexMatrix& data() const {return wigner_;}

        //initialization
        /**
         * @brief Initial state of the WDF - default is gaussian
         * @param x0,p0 - central part of the Gaussian
         * @param sigma_x, sigma_p - uncertainties in the Gaussian
         * @brief we can also initialize from function - initializeFromFunction()
         */
        void initializeGaussian(double x0, double p0, double sigma_x, double sigma_p);
        void initializeCoherentState(double x0, double p0, double sigma);
        void initializeFromFunction(const std::function<Complex(double, double)>& func);

        // expectation values
        double norm() const;
        double totalProbability() const;
        std::pair<double, double> expectedPosition() const;
        std::pair<double, double> expectedMomentum() const;

        //saving
        void saveToFile(const std::string& filename) const;

        // visualization data
        RealMatrix real() const;
        RealMatrix magnitude() const;

    private:
        void normalize();
};