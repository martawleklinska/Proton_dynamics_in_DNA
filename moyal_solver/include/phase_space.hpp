#pragma once
#include "moyal_setup.hpp"
#include<fftw3.h>

/**
 * @brief Phase space grid and fft
 */

 class PhaseSpace{
    private:
        MoyalConfig config_;
        std::vector<double> x_vec_;      //!< position grid
        std::vector<double> p_vec_;      //!< momentum grid
        std::vector<double> kx_vec_;     //!< position frequencies
        std::vector<double> kp_vec_;     //!< momentum frequencies

        RealMatrix X_, P_;              //!< Meshgrids
        RealMatrix KX_, KP_;            //!< Frequency meshgrids

        double dx_, dp_;                //!< grid spacing

        //fftw plans
        fftw_plan fft_x_forward_, fft_x_backward_;
        fftw_plan fft_p_forward_, fft_p_backward_;

        fftw_complex* work_x_;
        fftw_complex* work_p_;

    public:
        explicit PhaseSpace(const MoyalConfig& config);
        ~PhaseSpace();

        //grid access
        const std::vector<double>& x_grid() const { return x_vec_; }
        const std::vector<double>& p_grid() const { return p_vec_; }
        const RealMatrix& X() const { return X_ ;}
        const RealMatrix& P() const { return P_ ;}
        const RealMatrix& KX() const { return KX_ ;}
        const RealMatrix& KP() const { return KP_ ;}
        
        double dx() const {return dx_;}
        double dp() const {return dp_;}
        int gridX() const {return config_.gridX;}
        int gridP() const {return config_.gridP;}

        /**
         *  @brief fft operations: copy data to work array, execute FFT normalize for backward and copy back
        */ 
        void fft_x(ComplexMatrix& data, bool forward = true);
        void fft_p(ComplexMatrix& data, bool forward = true);
        /**
         * 
         */
        void fftshift_x(ComplexMatrix& data);
        void fftshift_p(ComplexMatrix& data);

    private:
        void initializeGrids();     //!< initialize: resize position, momentum and ffts grids
        void setupFFT();            //!< allocate momory for each plan
        void createMeshgrids();     //!< setup meshgrids - reside matrices X_, P_, KX_, KP_
};
