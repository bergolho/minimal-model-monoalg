#include <stdio.h>
#include "fhn_mod.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using modified FHN 1961 CPU model\n");

    uint32_t num_cells = solver->original_num_cells;
    solver->sv = (real*)malloc(NEQ*num_cells*sizeof(real));

    OMP(parallel for)
    for(uint32_t i = 0; i < num_cells; i++) {
        real *sv = &solver->sv[i * NEQ];

        sv[0] = 0.000000f; // Vm millivolt
        sv[1] = 0.000000f; // v dimensionless
    }
}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    uint32_t sv_id;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    OMP(parallel for private(sv_id))
    for (uint32_t i = 0; i < num_cells_to_solve; i++) {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = (uint32_t )i;

        for (int j = 0; j < num_steps; ++j) {
            solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i]);

        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  {

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current);

    for(int i = 0; i < NEQ; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current) {

    //State variables
    const real v = sv[0];
    const real w = sv[1];

    const real G = 1.5f;			// omega^-1 * cm^-2
    const real eta1 = 4.4f;			// dimensionless
    const real eta2 = 0.012f;		// dimensionless
    const real eta3 = 1.0f;			// dimensionless
    const real vth = 13.0f;			// mV
    const real vp = 100.0f;			// mV
    
    const real chi_ = 1.0e3;		// cm^-1
    const real Cm_ = 1.0e-3;		// mF * cm^-2


    rDY_[0] = ((1.0f / (Cm_ * chi_)) * ((-G * v * (1.0f - (v / vth)) * (1.0f - (v / vp))) + (-eta1 * v * w))) + stim_current;
    rDY_[1] = (eta2 * ((v / vp) - (eta3 * w)));


}

