#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "minimal_model.h"
//#include "functions.inc.c"
#include <stdio.h>

GET_CELL_MODEL_DATA(init_cell_model_data) {

   assert(cell_model);

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;

}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using Minimal Model Mixed CPU model\n");

// Get the mapping array
    uint32_t *mapping = NULL;

    if(solver->ode_extra_data)
    {
        mapping = (uint32_t*)solver->ode_extra_data;
    }
    else
    {
        log_error_and_exit("You need to specify a mask function when using a mixed model!\n");
    }

    uint32_t num_cells = solver->original_num_cells;

    solver->sv = (real*)malloc(NEQ*num_cells*sizeof(real));

    OMP(parallel for)
        for(uint32_t i = 0; i < num_cells; i++) {
            real *sv = &solver->sv[i * NEQ];
            sv[0] = 0.0f; //u
            sv[1] = 1.0f; //v
            sv[2] = 1.0f; //w
            sv[3] = 0.0f; //s
        }
}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    // Get the mapping array
    uint32_t *mapping = NULL;
    if(ode_solver->ode_extra_data)
    {
        mapping = (uint32_t*)ode_solver->ode_extra_data;
    }
    else
    {
        log_error_and_exit("You need to specify a mask function when using a mixed model!\n");
    }

    uint32_t sv_id;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    int i;

    OMP(parallel for private(sv_id))
    for (i = 0; i < num_cells_to_solve; i++) {
        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        for (int j = 0; j < num_steps; ++j)
        {
            // nao tem necessidade desses ifs, ja que passa o mapping[i] em vez de valores 0, 1 e 2.
            //if (mapping[i] == 0)          // ENDO
            //    solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], mapping[i]);
            //else if (mapping[i] == 1)     // MYO
            //    solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], mapping[i]);
            //else                          // EPI
            //    solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], mapping[i]);
            solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], mapping[i]);
        }
    }

}


void solve_model_ode_cpu(real dt, real *sv, real stim_current, int type_cell)  {

    assert(sv);

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt, type_cell);

    //THIS MODEL USES THE Rush Larsen Method TO SOLVE THE EDOS
    sv[0] = dt*rDY[0] + rY[0];
    sv[1] = dt*rDY[1] + rY[1];
    sv[2] = dt*rDY[2] + rY[2];
    sv[3] = dt*rDY[3] + rY[3];
}

// lado direito
void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, int type_cell) {

    //THIS IS THE STATE VECTOR THAT WE NEED TO SAVE IN THE STEADY STATE
    const real s_u = sv[0];
    const real s_v = sv[1];
    const real s_w = sv[2];
    const real s_s = sv[3];

    //////////////////////////////////////////

    const real u_o = 0.0;
    const real theta_v = 0.3;
    const real theta_w = 0.13;
    const real tau_vplus = 1.4506;
    const real tau_s1 = 2.7342;
    const real k_s = 2.0994;
    const real u_s = 0.9087;

    if (type_cell == 0) {        // ENDO
      //const real u_o = 0.0;
      const real u_u = 1.56;
      //const real theta_v = 0.3;
      //const real theta_w = 0.13;
      const real theta_vminus = 0.2;
      const real theta_o = 0.006;
      const real tau_v1minus = 75.0;
      const real tau_v2minus = 10.0;
      //const real tau_vplus = 1.4506;
      const real tau_w1minus = 6.0;
      const real tau_w2minus = 140.0;
      const real k_wminus = 200.0;
      const real u_wminus = 0.016;
      const real tau_wplus = 280.0;
      const real tau_fi = 0.1;
      const real tau_o1 = 470.0;
      const real tau_o2 = 6.0;
      const real tau_so1 = 40.0;
      const real tau_so2 = 1.2;
      const real k_so = 2.0;
      const real u_so = 0.65;
      //const real tau_s1 = 2.7342;
      const real tau_s2 = 2.0;
      //const real k_s = 2.0994;
      //const real u_s = 0.9087;
      const real tau_si = 2.9013;
      const real tau_winf = 0.0273;
      const real w_infstar = 0.78;
    }
    else if (type_cell == 1) {   // MYO
      //const real u_o = 0.0;
      const real u_u = 1.61;
      //const real theta_v = 0.3;
      //const real theta_w = 0.13;
      const real theta_vminus = 0.1;
      const real theta_o = 0.005;
      const real tau_v1minus = 80.0;
      const real tau_v2minus = 1.4506;
      //const real tau_vplus = 1.4506;
      const real tau_w1minus = 70.0;
      const real tau_w2minus = 8.0;
      const real k_wminus = 200.0;
      const real u_wminus = 0.016;
      const real tau_wplus = 280.0;
      const real tau_fi = 0.078;
      const real tau_o1 = 410.0;
      const real tau_o2 = 7.0;
      const real tau_so1 = 91.0;
      const real tau_so2 = 0.8;
      const real k_so = 2.1;
      const real u_so = 0.6;
      //const real tau_s1 = 2.7342;
      const real tau_s2 = 4.0;
      //const real k_s = 2.0994;
      //const real u_s = 0.9087;
      const real tau_si = 3.3849;
      const real tau_winf = 0.01;
      const real w_infstar = 0.5;
    }
    else {                       // EPI
      //const real u_o = 0.0;
      const real u_u = 1.55;
      //const real theta_v = 0.3;
      //const real theta_w = 0.13;
      const real theta_vminus = 0.006;
      const real theta_o = 0.006;
      const real tau_v1minus = 60.0;
      const real tau_v2minus = 1150.0;
      //const real tau_vplus = 1.4506;
      const real tau_w1minus = 60.0;
      const real tau_w2minus = 15.0;
      const real k_wminus = 65.0;
      const real u_wminus = 0.03;
      const real tau_wplus = 200.0;
      const real tau_fi = 0.11;
      const real tau_o1 = 400.0;
      const real tau_o2 = 6.0;
      const real tau_so1 = 30.0181;
      const real tau_so2 = 0.9957;
      const real k_so = 2.0458;
      const real u_so = 0.65;
      //const real tau_s1 = 2.7342;
      const real tau_s2 = 16.0;
      //const real k_s = 2.0994;
      //const real u_s = 0.9087;
      const real tau_si = 1.8875;
      const real tau_winf = 0.07;
      const real w_infstar = 0.94;
    }

        ustep = u;
        vstep = v;
        wstep = w;
        sstep = s;

        // Get du_dt, dv_dt, dw_dt and ds_dt
        //du_dt = - reaction_u(ustep, vstep, wstep, sstep) + I_stim;
        //
        const real H = 0.0;
        const real h_o = 0.0;
        const real h_w = 0.0;
        const real h_v_minus = 0.0;
        const real tau_vminus = 0.0;
        const real J_fi = 0.0;
        const real J_so = 0.0;
        const real J_si = 0.0;
        const real tau_o = 0.0;
        const real tau_so = 0.0;
        const real tau_s = 0.0;
        const real du_dt = 0.0;
        const real dv_dt = 0.0;
        const real dw_dt = 0.0;
        const real ds_dt = 0.0;
        const real v_inf = 0.0;
        const real w_inf = 0.0;
        const real tau_wminus = 0.0;

        if (u-theta_v > 0)
          H = 1.0;
        J_fi = -v * H * (u-theta_v) * (u_u-u) / tau_fi;
        if (u-theta_o > 0)
          h_o = 1.0;
        tau_o = (1.0 - h_o) * tau_o1 + (h_o * tau_o2);
        if (u-theta_w > 0)
          h_w = 1.0;
        tau_so = tau_so1 + (((tau_so2 - tau_so1) * (1.0 + tanh(k_so*(u - u_so)))) * 0.5);
        J_so = ((u-u_o) * (1.0 - h_w) / tau_o) + (h_w / tau_so);
        J_si = - h_w * w * s / tau_si;
        du_dt = - (J_fi + J_so + J_si) + I_stim;
        //

        //dv_dt = dvdt(ustep, vstep);
        //
        if (u < theta_vminus)
        {
            v_inf = 1.0;
        }
        if (u-theta_vminus > 0)
          h_v_minus = 1.0;
        tau_vminus = (1.0 - h_v_minus) * tau_v1minus + (h_v_minus * tau_v2minus);
        dv_dt = (1.0 - H) * (v_inf - v) / tau_vminus - (H * v / tau_vplus);
        //

        //dw_dt = dwdt(ustep, wstep);
        //
        w_inf = (1.0 - h_o) * (1.0 - (u/tau_winf)) + (h_o * w_infstar);
        tau_wminus = tau_w1minus + (((tau_w2minus - tau_w1minus) * (1.0 + tanh(k_wminus*(u - u_wminus)))) * 0.5);
        dw_dt = (1.0 - h_w) * (w_inf - w) / tau_wminus - (h_w * w / tau_wplus)
        //

        //ds_dt = dsdt(ustep, sstep);
        //
        tau_s = (1.0 - h_w) * tau_s1 + (h_w * tau_s2);
        ds_dt = (((1.0 + tanh(k_s*(u - u_s))) * 0.5) - s) / tau_s;
        //

        // Update u, v, w and s
        u = ustep + dt * du_dt;
        v = vstep + dt * dv_dt;
        w = wstep + dt * dw_dt;
        s = sstep + dt * ds_dt;

    //////////////////////////////////////////
}