//
//  main.cpp
//  MaxwellStefanDynamic
//
//  Created by Derek Harrison on 09/03/2022.
//

#include <iostream>

#include <stdio.h>
#include <cstdlib>
#include <time.h>

typedef struct set_up_params {
    double V;
    double d;
    double len;
} su_params_t;

typedef struct time_params {
    double to;
    double tf;
    double dt;
} t_params_t;

typedef struct physical_params {
    double ** D;
    double ct;
} p_params_t;

typedef struct bulb_composition {
    double * x_b1;
    double * x_b2;
} b_comp_t;

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

int main(int argc, char* argv[]) {

    // Number of components and grid resolution
    int num_components = 3;
    int num_grid_cells = 6;
    
    // Input concentration
    p_params_t physical_params;
    physical_params.ct = 1.0;
    
    // Diffusivities
    double D12 = 8.33e-5 * 3600; // units are (m2 / h)
    double D13 = 6.8e-5 * 3600; // units are (m2 / h)
    double D23 = 1.68e-5 * 3600; // units are (m2 / h)
    
    // Experimental setup parameters
    su_params_t set_up_params;
    set_up_params.V = 5e-4;
    set_up_params.d = 2e-3;
    set_up_params.len = 1e-2;
    
    // Time parameters
    t_params_t time_params;
    time_params.to = 0.0;
    time_params.tf = 20.0;
    time_params.dt = 1e-6; // This value has to be low for accuracy
    
    // Initial composition bulb1
    double xb10 = 0.501; // H2 fraction
    double xb11 = 0.499; // N2 fraction
    double xb12 = 1 - xb10 - xb11; // CO2 fraction
    
    // Initial composition bulb2
    double xb20 = 0; // H2 fraction
    double xb21 = 0.501; // N2 fraction
    double xb22 = 1 - xb20 - xb21; // CO2 fraction
    
    // Set up input
    physical_params.D = mat2D(num_components, num_components);
    
    double *** A = mat3D(num_grid_cells + 1, num_components, num_components);
    double *** L = mat3D(num_grid_cells + 1, num_components, num_components);
    double *** U = mat3D(num_grid_cells + 1, num_components, num_components);
    
    b_comp_t bulb_compositions;
    bulb_compositions.x_b1 = new double[num_components];
    bulb_compositions.x_b2 = new double[num_components];
    
    double ** tube_fracs = mat2D(num_grid_cells, num_components);
    double ** grad_vec = mat2D(num_grid_cells + 1, num_components);
    double ** flux_vec = mat2D(num_grid_cells + 1, num_components);
    double ** z_vec = mat2D(num_grid_cells + 1, num_components);
   
    bulb_compositions.x_b1[0] = xb10;
    bulb_compositions.x_b1[1] = xb11;
    bulb_compositions.x_b1[2] = xb12;
    
    bulb_compositions.x_b2[0] = xb20;
    bulb_compositions.x_b2[1] = xb21;
    bulb_compositions.x_b2[2] = xb22;
    
    physical_params.D[0][1] = D12;
    physical_params.D[0][2] = D13;
    physical_params.D[1][2] = D23;
    
    physical_params.D[1][0] = D12;
    physical_params.D[2][0] = D13;
    physical_params.D[2][1] = D23;
    
    // Initialize tube region
    init_tube(num_grid_cells, num_components, tube_fracs);
    
    // Perform simulation
    double t = time_params.to;
    
    while(t < time_params.tf) {

        compute_linear_system(grad_vec,
                              tube_fracs,
                              num_grid_cells,
                              num_components,
                              bulb_compositions,
                              physical_params,
                              set_up_params,
                              A);
        
        for(int g = 0; g < num_grid_cells + 1; ++g) {
            lu_decomposition(A[g], num_components, L[g], U[g]);
            
            solve_Lz(L[g], num_components, grad_vec[g], z_vec[g]);

            solve_Uy(U[g], num_components, z_vec[g], flux_vec[g]);
        }
        
        // Compute flux of component n
        compute_flux_n(num_grid_cells, flux_vec, num_components);
        
        update_composition(num_grid_cells,
                           num_components,
                           flux_vec,
                           time_params,
                           physical_params,
                           bulb_compositions,
                           set_up_params,
                           tube_fracs);
        
        t = t + time_params.dt;
    }
    
    // Print results
    std::cout << "bulb 1" << std::endl;
    for(int c = 0; c < num_components; ++c)
        std::cout << "Fraction of component " << c << ": " << bulb_compositions.x_b1[c] << std::endl;
    
    std::cout << "bulb 2" << std::endl;
    for(int c = 0; c < num_components; ++c)
        std::cout << "Fraction of component " << c << ": " << bulb_compositions.x_b2[c] << std::endl;
    
    return 0;
}

