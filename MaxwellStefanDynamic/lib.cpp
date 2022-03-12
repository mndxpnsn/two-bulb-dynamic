//
//  lib.cpp
//  MaxwellStefanDynamic
//
//  Created by mndx on 10/03/2022.
//

#include <iostream>

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

void free_mat2D(double ** mat, int nx) {
    
    for(int i = 0; i < nx; ++i) {
        delete [] mat[i];
    }
    
    delete [] mat;
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

void free_mat3D(double *** mat, int nx, int ny) {
    
    for(int i = 0; i < nx; ++i) {
        for(int j = 0; j < ny; ++j) {
            delete [] mat[i][j];
        }
        
        delete [] mat[i];
    }
    
    delete [] mat;
}

void lu_decomposition(double ** A, int N, double ** L, double ** U) {

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

void compute_linear_system(int n,
                           int ng,
                           s_data_t sim_data,
                           e_data_t experiment_data,
                           b_comp_t b_comp,
                           double *** a) {
    
    double dz = experiment_data.su_params.len / ng;
    double ct = experiment_data.p_params.ct;
    
    // Left most node
    for(int c = 0; c < n; ++c) {
        sim_data.grad_vec[0][c] = -ct * (sim_data.tube_fracs[0][c] - b_comp.x_b1[c]) / (0.5 * dz);
    }
    
    for(int i = 0; i < n; ++i) {
        double aii = 0.0;
        double diff = 1.0;
        for(int k = 0; k < n; ++k) {
            diff = diff - b_comp.x_b1[k];
            if(i != k) {
                double xk_loc = b_comp.x_b1[k];
                aii = aii + xk_loc / experiment_data.p_params.D[i][k];
            }
        }
        aii = aii + b_comp.x_b1[i] / experiment_data.p_params.D[i][n - 1];
        a[0][i][i] = aii;
    }
    
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            if(i != j) {
                a[0][i][j] = b_comp.x_b1[i] *
                            (1.0 / experiment_data.p_params.D[0][n - 1] - 1.0 / experiment_data.p_params.D[i][j]);
            }
        }
    }
    
    // Middle nodes
    for(int g = 1; g < ng; ++g) {
        for(int c = 0; c < n; ++c) {
            sim_data.grad_vec[g][c] = -ct * (sim_data.tube_fracs[g][c] - sim_data.tube_fracs[g - 1][c]) / dz;
        }
    }
    
    for(int g = 1; g < ng; ++g) {
        for(int i = 0; i < n; ++i) {
            double aii = 0.0;
            double diff = 1.0;
            for(int k = 0; k < n; ++k) {
                diff = diff - sim_data.tube_fracs[g][k];
                if(i != k) {
                    double xk_loc = (sim_data.tube_fracs[g][k] + sim_data.tube_fracs[g - 1][k]) / 2;
                    aii = aii + xk_loc / experiment_data.p_params.D[i][k];
                }
            }
            aii = aii + sim_data.tube_fracs[g][i] / experiment_data.p_params.D[i][n - 1];
            a[g][i][i] = aii;
        }
        
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                if(i != j) {
                    a[g][i][j] = sim_data.tube_fracs[g][i] *
                                (1.0 / experiment_data.p_params.D[i][n - 1] - 1.0 / experiment_data.p_params.D[i][j]);
                }
            }
        }
    }
    
    // Node n
    for(int c = 0; c < n; ++c) {
        sim_data.grad_vec[ng][c] = -ct * (b_comp.x_b2[c] - sim_data.tube_fracs[ng - 1][c]) / (0.5 * dz);
    }
    
    for(int i = 0; i < n; ++i) {
        double aii = 0.0;
        double diff = 1.0;
        for(int k = 0; k < n; ++k) {
            diff = diff - b_comp.x_b2[k];
            if(i != k) {
                double xk_loc = b_comp.x_b2[k];
                aii = aii + xk_loc / experiment_data.p_params.D[i][k];
            }
        }
        aii = aii + b_comp.x_b2[i] / experiment_data.p_params.D[i][n - 1];
        a[ng][i][i] = aii;
    }
    
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            if(i != j) {
                a[ng][i][j] = b_comp.x_b2[i] *
                            (1.0 / experiment_data.p_params.D[i][n - 1] - 1.0 / experiment_data.p_params.D[i][j]);
            }
        }
    }
}

void update_composition(int ng,
                        int n,
                        s_data_t sim_data,
                        e_data_t experiment_data,
                        b_comp_t b_comps) {
    
    double dt = experiment_data.time_params.dt;
    double d = experiment_data.su_params.d;
    double V = experiment_data.su_params.V;
    double ct = experiment_data.p_params.ct;
    double dz = experiment_data.su_params.len / ng;
    double Ac = 3.14 * d * d / 4.0;
    
    // Bulb1
    double diff = 1.0;
    for(int c = 0; c < n - 1; ++c) {
        b_comps.x_b1[c] = b_comps.x_b1[c] - sim_data.flux_vec[0][c] * dt * Ac / (V * ct);
        diff = diff - b_comps.x_b1[c];
    }
    
    // Update mole fraction of component n in bulb 1
    b_comps.x_b1[n - 1] = diff;
    
    // Tube
    for(int g = 0; g < ng; ++g) {
        diff = 1.0;
        for(int c = 0; c < n - 1; ++c) {
            sim_data.tube_fracs[g][c] = sim_data.tube_fracs[g][c] -
            (sim_data.flux_vec[g + 1][c] - sim_data.flux_vec[g][c]) * dt / (ct * dz);
            
            diff = diff - sim_data.tube_fracs[g][c];
        }
        
        //Update mole fraction of component n in the tube
        sim_data.tube_fracs[g][n - 1] = diff;
    }
    
    // Bulb2
    diff = 1.0;
    for(int c = 0; c < n - 1; ++c) {
        b_comps.x_b2[c] = b_comps.x_b2[c] + sim_data.flux_vec[ng][c] * dt * Ac / (V * ct);
        diff = diff - b_comps.x_b2[c];
    }
    
    // Update mole fraction of component n in bulb 2
    b_comps.x_b2[n - 1] = diff;
}

void compute_flux_n(int n, int ng, s_data_t sim_data) {
    //Compute flux of component n
    for(int g = 0; g < ng; ++g) {
        double Jn = 0.0;
        for(int c = 0; c < n - 1; ++c) {
            Jn = Jn - sim_data.flux_vec[g][c];
        }
        sim_data.flux_vec[g][n - 1] = Jn;
    }
}

void init_tube(int ng, int n, double ** x) {
    for(int g = 0; g < ng; ++g) {
        for(int c = 0; c < n; ++c) {
            x[g][c] = 1.0 / n;
        }
    }
}

void compute_compositions(int num_components,
                          int num_grid_cells,
                          e_data_t experiment_data,
                          b_comp_t bulb_compositions) {
    
    // LU decomposition data
    double *** A = mat3D(num_grid_cells + 1, num_components, num_components);
    double *** L = mat3D(num_grid_cells + 1, num_components, num_components);
    double *** U = mat3D(num_grid_cells + 1, num_components, num_components);
    
    // Simulatiton data
    s_data_t sim_data;
    sim_data.tube_fracs = mat2D(num_grid_cells, num_components);
    sim_data.grad_vec = mat2D(num_grid_cells + 1, num_components);
    sim_data.flux_vec = mat2D(num_grid_cells + 1, num_components);
    sim_data.z_vec = mat2D(num_grid_cells + 1, num_components);
    
    // Initialize tube region
    init_tube(num_grid_cells, num_components, sim_data.tube_fracs);
    
    // Perform simulation
    double t = experiment_data.time_params.to;
    double tf = experiment_data.time_params.tf;
    double dt = experiment_data.time_params.dt;
    
    while(t < tf) {

        compute_linear_system(num_components,
                              num_grid_cells,
                              sim_data,
                              experiment_data,
                              bulb_compositions,
                              A);
        
        for(int g = 0; g < num_grid_cells + 1; ++g) {
            
            lu_decomposition(A[g],
                             num_components,
                             L[g],
                             U[g]);
            
            solve_Lz(L[g],
                     num_components,
                     sim_data.grad_vec[g],
                     sim_data.z_vec[g]);

            solve_Uy(U[g],
                     num_components,
                     sim_data.z_vec[g],
                     sim_data.flux_vec[g]);
        }
        
        // Compute flux of component n
        compute_flux_n(num_components,
                       num_grid_cells,
                       sim_data);
        
        update_composition(num_grid_cells,
                           num_components,
                           sim_data,
                           experiment_data,
                           bulb_compositions);
        
        t = t + dt;
    }
    
    // Deallocate data
    free_mat3D(A, num_grid_cells + 1, num_components);
    free_mat3D(U, num_grid_cells + 1, num_components);
    free_mat3D(L, num_grid_cells + 1, num_components);
    
    free_mat2D(sim_data.tube_fracs, num_grid_cells);
    free_mat2D(sim_data.grad_vec, num_grid_cells + 1);
    free_mat2D(sim_data.flux_vec, num_grid_cells + 1);
    free_mat2D(sim_data.z_vec, num_grid_cells + 1);
}

void print_results(int num_components, b_comp_t bulb_compositions) {
    
    std::cout << "bulb 1" << std::endl;
    for(int c = 0; c < num_components; ++c) {
        std::cout << "Fraction of component " << c << ": ";
        std::cout << bulb_compositions.x_b1[c] << std::endl;
    }
    
    std::cout << "bulb 2" << std::endl;
    for(int c = 0; c < num_components; ++c) {
        std::cout << "Fraction of component " << c << ": ";
        std::cout << bulb_compositions.x_b2[c] << std::endl;
    }
}
