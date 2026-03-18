#include "phase_space.hpp"
#include <cmath>
#include <stdexcept>

const int FFTW_N_THREADS = 8;

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
    int res = fftw_init_threads();
    if (res == 0){
        throw std::runtime_error("FFTW thread init failed");
    }
    work_x_ = fftw_alloc_complex(config_.gridX);
    work_p_ = fftw_alloc_complex(config_.gridP);
    
    if (!work_x_ || !work_p_) {
        throw std::runtime_error("Failed to allocate FFTW memory");
    }
    fftw_plan_with_nthreads(FFTW_N_THREADS);
    fft_x_forward_ = fftw_plan_dft_1d(config_.gridX, work_x_, work_x_, 
                                     FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan_with_nthreads(FFTW_N_THREADS);
    fft_x_backward_ = fftw_plan_dft_1d(config_.gridX, work_x_, work_x_, 
                                      FFTW_BACKWARD, FFTW_MEASURE);
    
    fftw_plan_with_nthreads(FFTW_N_THREADS);
    fft_p_forward_ = fftw_plan_dft_1d(config_.gridP, work_p_, work_p_, 
                                     FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan_with_nthreads(FFTW_N_THREADS);
    fft_p_backward_ = fftw_plan_dft_1d(config_.gridP, work_p_, work_p_, 
                                      FFTW_BACKWARD, FFTW_MEASURE);
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

void PhaseSpace::applyMomentumMask(ComplexMatrix& data, double mask_fraction) {
    int nx = config_.gridX;
    int np = config_.gridP;
    
    // Szerokość strefy tłumienia (np. 10% siatki z każdej strony)
    int mask_width = static_cast<int>(mask_fraction * np);
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            double factor = 1.0;
            
            // Lewa strona (j bliskie 0)
            if (j < mask_width) {
                double t = static_cast<double>(j) / mask_width;
                factor = 0.5 * (1.0 - std::cos(M_PI * t)); // cosine window
            }
            // Prawa strona (j bliskie np-1)
            else if (j >= np - mask_width) {
                double t = static_cast<double>(np - 1 - j) / mask_width;
                factor = 0.5 * (1.0 - std::cos(M_PI * t));
            }
            
            data[i][j] *= factor;
        }
    }
}

void PhaseSpace::applyPositionMask(ComplexMatrix& data, double mask_fraction) {
    int nx = config_.gridX;
    int np = config_.gridP;
    int mask_width = static_cast<int>(mask_fraction * nx);
    
    for (int i = 0; i < nx; ++i) {
        double factor = 1.0;
        if (i < mask_width) {
            double t = static_cast<double>(i) / mask_width;
            factor = 0.5 * (1.0 - std::cos(M_PI * t));
        } else if (i >= nx - mask_width) {
            double t = static_cast<double>(nx - 1 - i) / mask_width;
            factor = 0.5 * (1.0 - std::cos(M_PI * t));
        }
        for (int j = 0; j < np; ++j) {
            data[i][j] *= factor;
        }
    }
}