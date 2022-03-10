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

#include "lib.hpp"
#include "user_types.h"

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

