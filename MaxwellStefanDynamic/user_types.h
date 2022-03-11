//
//  user_types.h
//  MaxwellStefanDynamic
//
//  Created by mndx on 10/03/2022.
//

#ifndef user_types_h
#define user_types_h

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

#endif /* user_types_h */
