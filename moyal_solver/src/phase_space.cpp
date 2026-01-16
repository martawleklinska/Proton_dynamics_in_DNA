#include "phase_space.hpp"
#include <cmath>
#include <stdexcept>

PhaseSpace::PhaseSpace(const MoyalConfig& config) : config_(config) {
    initializeGrids();
    createMeshgrids();
    setupFFT();
}

PhaseSpace::~PhaseSpace() {
    if (fft_x_forward_) fftw_destroy_plan(fft_x_forward_);
    if (fft_x_backward_) fftw_destroy_plan(fft_x_backward_);
    if (fft_p_forward_) fftw_destroy_plan(fft_p_forward_);
    if (fft_p_backward_) fftw_destroy_plan(fft_p_backward_);
    if (work_x_) fftw_free(work_x_);
    if (work_p_) fftw_free(work_p_);
}

void PhaseSpace::initializeGrids() {
    x_vec_.resize(config_.gridX);
    dx_ = config_.ampX / config_.gridX;
    for (int i = 0; i < config_.gridX; ++i) {
        x_vec_[i] = -config_.ampX/2.0 + i * dx_;
    }
    
    p_vec_.resize(config_.gridP);
    dp_ = config_.ampP / config_.gridP;
    for (int i = 0; i < config_.gridP; ++i) {
        p_vec_[i] = -config_.ampP/2.0 + i * dp_;
    }

    kx_vec_.resize(config_.gridX);
    kp_vec_.resize(config_.gridP);
    
    for (int i = 0; i < config_.gridX; ++i) {
        int n = (i < config_.gridX/2) ? i : i - config_.gridX;
        kx_vec_[i] = 2.0 * M_PI * n / (config_.gridX * dx_);
    }
    
    for (int i = 0; i < config_.gridP; ++i) {
        int n = (i < config_.gridP/2) ? i : i - config_.gridP;
        kp_vec_[i] = 2.0 * M_PI * n / (config_.gridP * dp_);
    }
}

void PhaseSpace::createMeshgrids() {
    X_.resize(config_.gridX, std::vector<double>(config_.gridP));
    P_.resize(config_.gridX, std::vector<double>(config_.gridP));
    KX_.resize(config_.gridX, std::vector<double>(config_.gridP));
    KP_.resize(config_.gridX, std::vector<double>(config_.gridP));
    
    for (int i = 0; i < config_.gridX; ++i) {
        for (int j = 0; j < config_.gridP; ++j) {
            X_[i][j] = x_vec_[i];
            P_[i][j] = p_vec_[j];
            KX_[i][j] = kx_vec_[i];
            KP_[i][j] = kp_vec_[j];
        }
    }
}

void PhaseSpace::setupFFT() {
    work_x_ = fftw_alloc_complex(config_.gridX);
    work_p_ = fftw_alloc_complex(config_.gridP);
    
    if (!work_x_ || !work_p_) {
        throw std::runtime_error("Failed to allocate FFTW memory");
    }
    
    fft_x_forward_ = fftw_plan_dft_1d(config_.gridX, work_x_, work_x_, 
                                     FFTW_FORWARD, FFTW_ESTIMATE);
    fft_x_backward_ = fftw_plan_dft_1d(config_.gridX, work_x_, work_x_, 
                                      FFTW_BACKWARD, FFTW_ESTIMATE);
    fft_p_forward_ = fftw_plan_dft_1d(config_.gridP, work_p_, work_p_, 
                                     FFTW_FORWARD, FFTW_ESTIMATE);
    fft_p_backward_ = fftw_plan_dft_1d(config_.gridP, work_p_, work_p_, 
                                      FFTW_BACKWARD, FFTW_ESTIMATE);
}

void PhaseSpace::fft_x(ComplexMatrix& data, bool forward) {
    for (int j = 0; j < config_.gridP; ++j) {
        for (int i = 0; i < config_.gridX; ++i) {
            work_x_[i][0] = data[i][j].real();
            work_x_[i][1] = data[i][j].imag();
        }
        
        if (forward) {
            fftw_execute(fft_x_forward_);
        } else {
            fftw_execute(fft_x_backward_);
            for (int i = 0; i < config_.gridX; ++i) {
                work_x_[i][0] /= config_.gridX;
                work_x_[i][1] /= config_.gridX;
            }
        }
        
        for (int i = 0; i < config_.gridX; ++i) {
            data[i][j] = Complex(work_x_[i][0], work_x_[i][1]);
        }
    }
}

void PhaseSpace::fft_p(ComplexMatrix& data, bool forward) {
    for (int i = 0; i < config_.gridX; ++i) {
        for (int j = 0; j < config_.gridP; ++j) {
            work_p_[j][0] = data[i][j].real();
            work_p_[j][1] = data[i][j].imag();
        }
        
        if (forward) {
            fftw_execute(fft_p_forward_);
        } else {
            fftw_execute(fft_p_backward_);
            for (int j = 0; j < config_.gridP; ++j) {
                work_p_[j][0] /= config_.gridP;
                work_p_[j][1] /= config_.gridP;
            }
        }
        
        for (int j = 0; j < config_.gridP; ++j) {
            data[i][j] = Complex(work_p_[j][0], work_p_[j][1]);
        }
    }
}

void PhaseSpace::fftshift_x(ComplexMatrix& data) {
    ComplexMatrix temp = data;
    int nx = config_.gridX;
    int np = config_.gridP;
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            int new_i = (i + nx/2) % nx;
            data[new_i][j] = temp[i][j];
        }
    }
}

void PhaseSpace::fftshift_p(ComplexMatrix& data) {
    ComplexMatrix temp = data;
    int nx = config_.gridX;
    int np = config_.gridP;
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            int new_j = (j + np/2) % np;
            data[i][new_j] = temp[i][j];
        }
    }
}
