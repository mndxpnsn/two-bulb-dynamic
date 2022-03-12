//
//  lib.hpp
//  MaxwellStefanDynamic
//
//  Created by mndx on 10/03/2022.
//

#ifndef lib_hpp
#define lib_hpp

#include <stdio.h>

#include "user_types.h"

double ** mat2D(int nx, int ny);

void free_mat2D(double ** mat, int nx);

double *** mat3D(int ng, int nx, int ny);

void lu_decomposition(double ** A, int N, double ** L, double ** U);

void solve_Uy(double ** U, int n, double * z, double * y);

void solve_Lz(double ** L, int n, double * x, double * z);

void mat_vec_mult(double ** A, int n, double * x, double * y);

void mat_mat_mult(double ** A, int n, double ** B, double ** C);

void compute_linear_system(double ** b_vec,
                           double ** x_vec,
                           int ng,
                           int n,
                           b_comp_t b_comp,
                           p_params_t p_params,
                           su_params_t su_params,
                           double *** a);

void update_composition(int ng,
                        int n,
                        double ** J,
                        t_params_t time_params,
                        p_params_t p_params,
                        b_comp_t b_comps,
                        su_params_t su_params,
                        double ** tube_fracs);

void compute_flux_n(int ng, double ** J, int n);

void init_tube(int ng, int n, double ** x);

void compute_compositions(int num_components,
                          int num_grid_cells,
                          e_data_t experiment_data,
                          b_comp_t bulb_compositions);

void print_results(int num_components, b_comp_t bulb_compositions);

#endif /* lib_hpp */
