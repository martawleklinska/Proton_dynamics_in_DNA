#pragma once
#include"moyal_setup.hpp"
#include"phase_space.hpp"
#include"propagators.hpp"
#include"wigner_df.hpp"

/**
 * @brief main Moyal equation solver using Strang splittling
 */

 class MoyalSolver{
    private:
        /**
         * @brief from earlier classes we take the config, phase space config, wdf and propagators
         */
        MoyalConfig config_;
        PhaseSpace phase_space_;
        WDF wigner_;

        std::unique_ptr<KineticPropagator> kinetic_prop_; 
        std::unique_ptr<PotentialPropagator> potential_prop_;

        int current_step_;
        double current_time_;

    public:
        explicit MoyalSolver(const MoyalConfig& config, std::unique_ptr<Potential> potential);

        // main evolution methods: initial condition from some function or in the general case - gaussian
        void setInitialCondition(const std::function<Complex(double, double)>& init_func);
        void setInitialGaussian(double x0, double p0, double sigma_x, double sigma_p);

        /**
         * @brief function evolve() saves current step and calles evolveOneStep()
         * evolveOneStep() applies strangSplittingStep()
         */
        void evolve(int steps = -1); // -1 means config_.timeSteps
        void evolveOneStep();
        // fft_p not with booleans

        //access classes functions in main
        const WDF& getWigner() const { return wigner_ ;}
        const PhaseSpace& getPhaseSpace() const { return phase_space_;}
        double getCurrentTime() const { return current_time_; }
        int getCurrentStep() const { return current_step_;}

        //analysis
        void computeExpectationValues(double& mean_x, double& mean_p,
                                        double& sigma_x, double& sigma_p) const;
        
        /**
         * @brief Calculate nonclassicality parameter \delta(t)
         * @return \delta(t) = \int |\varrho(x,p,t)| dx dp - 1
         */
        double calculateNonclassicalityParameter() const;

        // saving
        void saveState(const std::string& filename) const;

    private:
        /**
         * @brief Strang splitting: exp(-i*dt/2*T) exp(-i*dt*V) exp(-i*dt/2*T)
         * We take half kinetic step, full potential and then again hald kinetic 
         * -> step means transform to p take a step and transform back
         */
        void strangSplittingStep();
        void outputData();
 }; 

