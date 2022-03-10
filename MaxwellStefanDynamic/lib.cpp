//
//  lib.cpp
//  MaxwellStefanDynamic
//
//  Created by mndx on 10/03/2022.
//

#include "lib.hpp"
#include "user_types.h"

double ** mat2D(int nx, int ny) {

    double ** f = new double * [nx];

    for(int i = 0; i < nx; ++i) {
        f[i] = new double[ny];
        for(int j = 0; j < ny; ++j) {
            f[i][j] = 1e-13;
        }
    }

    return f;
}

double *** mat3D(int ng, int nx, int ny) {

    double *** f = new double ** [ng];

    for(int i = 0; i < ng; ++i) {
        f[i] = new double * [nx];
    }
    
    for(int g = 0; g < ng; ++g) {
        for(int i = 0; i < nx; ++i) {
            f[g][i] = new double[ny];
            for(int j = 0; j < ny; ++j) {
                f[g][i][j] = 1e-13;
            }
        }
        
    }

    return f;
}

void lu_decomposition(double ** A,
                      int N,
                      double ** L,
                      double ** U) {

    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            if(i == j) {
                L[i][j] = 1.0;
            }
            else {
                L[i][j] = 0.0;
            }
        }
    }

    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            if(i > j) {
                U[i][j] = 0.0;
            }
        }
    }

    for(int k = 0; k < N; ++k) {
        U[k][k] = A[k][k];
        for(int i = k + 1; i < N; ++i) {
            L[i][k] = A[i][k] / (U[k][k] + 1e-13);
            U[k][i] = A[k][i];
        }
        for(int i = k + 1; i < N; ++i) {
            for(int j = k + 1; j < N; ++j) {
                A[i][j] = A[i][j] - L[i][k]*U[k][j];
            }
        }
    }
}

void solve_Uy(double ** U, int n, double * z, double * y) {
    
    for(int i = n - 1; i >= 0; --i) {
        double yi = z[i];
        for(int j = i + 1; j < n; ++j) {
            yi = yi - U[i][j] * y[j];
        }
        yi = yi / (U[i][i] + 1e-13);
        y[i] = yi;
    }
}

void solve_Lz(double ** L, int n, double * x, double * z) {
    
    for(int i = 0; i < n; ++i) {
        double zi = x[i];
        for(int j = 0; j < i; ++j) {
            zi = zi - L[i][j] * z[j];
        }
        zi = zi / (L[i][i] + 1e-13);
        z[i] = zi;
    }
}

void mat_vec_mult(double ** A, int n, double * x, double * y) {
    
    for(int i = 0; i < n; ++i) {
        double y_elem = 0.0;
        for(int k = 0; k < n; ++k) {
            y_elem = y_elem + A[i][k] * x[k];
        }
        y[i] = y_elem;
    }
}

void mat_mat_mult(double ** A, int n, double ** B, double ** C) {
    
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            double elem = 0.0;
            for(int k = 0; k < n; ++k) {
                elem = elem + A[i][k] * B[k][j];
            }
            C[i][j] = elem;
        }
    }
}

void compute_linear_system(double ** b_vec,
                           double ** x_vec,
                           int ng,
                           int n,
                           b_comp_t b_comp,
                           p_params_t p_params,
                           su_params_t su_params,
                           double *** a) {
    
    double dz = su_params.len / ng;
    double ct = p_params.ct;
    
    //0
    for(int c = 0; c < n; ++c) {
        b_vec[0][c] = -ct * (x_vec[0][c] - b_comp.x_b1[c]) / (0.5 * dz);
    }
    
    for(int i = 0; i < n; ++i) {
        double aii = 0.0;
        double diff = 1.0;
        for(int k = 0; k < n; ++k) {
            diff = diff - b_comp.x_b1[k];
            if(i != k) {
                double xk_loc = b_comp.x_b1[k];
                aii = aii + xk_loc / p_params.D[i][k];
            }
        }
        aii = aii + b_comp.x_b1[i] / p_params.D[i][n - 1];
        a[0][i][i] = aii;
    }
    
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            if(i != j) {
                a[0][i][j] = b_comp.x_b1[i] * (1.0 / p_params.D[0][n - 1] - 1.0 / p_params.D[i][j]);
            }
        }
    }
    
    // mid
    for(int g = 1; g < ng; ++g) {
        for(int c = 0; c < n; ++c) {
            b_vec[g][c] = -ct * (x_vec[g][c] - x_vec[g - 1][c]) / dz;
        }
    }
    
    for(int g = 1; g < ng; ++g) {
        for(int i = 0; i < n; ++i) {
            double aii = 0.0;
            double diff = 1.0;
            for(int k = 0; k < n; ++k) {
                diff = diff - x_vec[g][k];
                if(i != k) {
                    double xk_loc = (x_vec[g][k] + x_vec[g - 1][k]) / 2;
                    aii = aii + xk_loc / p_params.D[i][k];
                }
            }
            aii = aii + x_vec[g][i] / p_params.D[i][n - 1];
            a[g][i][i] = aii;
        }
        
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                if(i != j) {
                    a[g][i][j] = x_vec[g][i] * (1.0 / p_params.D[i][n - 1] - 1.0 / p_params.D[i][j]);
                }
            }
        }
    }
    
    // n
    for(int c = 0; c < n; ++c) {
        b_vec[ng][c] = -ct * (b_comp.x_b2[c] - x_vec[ng - 1][c]) / (0.5 * dz);
    }
    
    for(int i = 0; i < n; ++i) {
        double aii = 0.0;
        double diff = 1.0;
        for(int k = 0; k < n; ++k) {
            diff = diff - b_comp.x_b2[k];
            if(i != k) {
                double xk_loc = b_comp.x_b2[k];
                aii = aii + xk_loc / p_params.D[i][k];
            }
        }
        aii = aii + b_comp.x_b2[i] / p_params.D[i][n - 1];
        a[ng][i][i] = aii;
    }
    
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            if(i != j) {
                a[ng][i][j] = b_comp.x_b2[i] * (1.0 / p_params.D[i][n - 1] - 1.0 / p_params.D[i][j]);
            }
        }
    }
}

void update_composition(int ng,
                        int n,
                        double ** J,
                        t_params_t time_params,
                        p_params_t p_params,
                        b_comp_t b_comps,
                        su_params_t su_params,
                        double ** tube_fracs) {
    
    double dt = time_params.dt;
    double d = su_params.d;
    double V = su_params.V;
    double ct = p_params.ct;
    double dz = su_params.len / ng;
    
    // Bulb1
    double diff = 1.0;
    for(int c = 0; c < n - 1; ++c) {
        b_comps.x_b1[c] = b_comps.x_b1[c] - J[0][c] * dt * 3.14 * d * d / 4.0 / (V * ct);
        diff = diff - b_comps.x_b1[c];
    }
    
    // Update mole fraction of component n in bulb 1
    b_comps.x_b1[n - 1] = diff;
    
    // Tube
    for(int g = 0; g < ng; ++g) {
        diff = 1.0;
        for(int c = 0; c < n - 1; ++c) {
            tube_fracs[g][c] = tube_fracs[g][c] - (J[g + 1][c] - J[g][c]) * dt / (ct * dz);
            diff = diff - tube_fracs[g][c];
        }
        
        //Update mole fraction of component n in the tube
        tube_fracs[g][n - 1] = diff;
    }
    
    // Bulb2
    diff = 1.0;
    for(int c = 0; c < n - 1; ++c) {
        b_comps.x_b2[c] = b_comps.x_b2[c] + J[ng][c] * dt * 3.14 * d * d / 4.0 / (V * ct);
        diff = diff - b_comps.x_b2[c];
    }
    
    // Update mole fraction of component n in bulb 2
    b_comps.x_b2[n - 1] = diff;
}

void compute_flux_n(int ng, double ** J, int n) {
    //Compute flux of component n
    for(int g = 0; g < ng; ++g) {
        double Jn = 0.0;
        for(int c = 0; c < n - 1; ++c) {
            Jn = Jn - J[g][c];
        }
        J[g][n - 1] = Jn;
    }
}

void init_tube(int ng, int n, double ** x) {
    for(int g = 0; g < ng; ++g) {
        for(int c = 0; c < n; ++c) {
            x[g][c] = 1.0 / n;
        }
    }
}
