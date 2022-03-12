//
//  main.cpp
//  MaxwellStefanDynamic
//
//  Created by mndx on 09/03/2022.
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
    
    // Total concentration
    p_params_t physical_params;
    physical_params.ct = 1.0; // Total concentration
    
    // Diffusivities
    double D12 = 8.33e-5 * 3600; // units are (m2 / h)
    double D13 = 6.8e-5 * 3600; // units are (m2 / h)
    double D23 = 1.68e-5 * 3600; // units are (m2 / h)
    
    // Experimental setup parameters
    su_params_t set_up_params;
    set_up_params.V = 5e-4; // Volume of the compartments (m3)
    set_up_params.d = 2e-3; // Diameter of tube (m)
    set_up_params.len = 1e-2; // Length of tube (m)
    
    // Time parameters
    t_params_t time_params;
    time_params.to = 0.0; // Initial time (h)
    time_params.tf = 20.0; // Final time (h)
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
    b_comp_t bulb_compositions;
    bulb_compositions.x_b1 = new double[num_components];
    bulb_compositions.x_b2 = new double[num_components];
    
    bulb_compositions.x_b1[0] = xb10;
    bulb_compositions.x_b1[1] = xb11;
    bulb_compositions.x_b1[2] = xb12;
    
    bulb_compositions.x_b2[0] = xb20;
    bulb_compositions.x_b2[1] = xb21;
    bulb_compositions.x_b2[2] = xb22;
    
    physical_params.D = mat2D(num_components, num_components);
    physical_params.D[0][1] = D12;
    physical_params.D[0][2] = D13;
    physical_params.D[1][2] = D23;
    
    physical_params.D[1][0] = D12;
    physical_params.D[2][0] = D13;
    physical_params.D[2][1] = D23;
    
    e_data_t experiment_data;
    experiment_data.p_params = physical_params;
    experiment_data.time_params = time_params;
    experiment_data.su_params = set_up_params;
    
    // Perform simulation of dynamic two-bulb experiment
    compute_compositions(num_components,
                         num_grid_cells,
                         experiment_data,
                         bulb_compositions);
    
    // Print results
    print_results(num_components, bulb_compositions);
    
    // Deallocate data
    delete [] bulb_compositions.x_b1;
    delete [] bulb_compositions.x_b2;
    free_mat2D(physical_params.D, num_components);
    
    return 0;
}

