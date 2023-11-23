#include "../../gpu_utils/gpu_utils.h"
#include "../../monodomain/constants.h"
#include <stddef.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "minimal_model.h"
#include <stdio.h>

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    log_info("Using Minimal Model Mixed GPU\n");

    uint32_t num_volumes = solver->original_num_cells;

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(real);

    check_cuda_error(cudaMallocPitch((void **) &(solver->sv), &pitch_h, size, (size_t )NEQ));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));

    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(solver->sv, NULL, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitch_h;
}

extern "C" SOLVE_MODEL_ODES(solve_model_odes_gpu) {

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    // execution configuration
    const int GRID  = ((int)num_cells_to_solve + BLOCK_SIZE - 1)/BLOCK_SIZE;


    size_t stim_currents_size = sizeof(real)*num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t)*num_cells_to_solve;

    real *stims_currents_device;
    check_cuda_error(cudaMalloc((void **) &stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));



    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **) &cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }

    // real *fibrosis_device;
    // real *fibs = NULL;
    // real fibs_size = num_cells_to_solve*sizeof(real);

    // bool deallocate = true;

    // check_cuda_error(cudaMalloc((void **) NULL, fibs_size));
    // check_cuda_error(cudaMemcpy(fibrosis_device, fibs, fibs_size, cudaMemcpyHostToDevice));

    uint32_t *mapping = NULL;
    uint32_t *mapping_device = NULL;
    if(ode_solver->ode_extra_data) 
    {
        mapping = (uint32_t*)ode_solver->ode_extra_data;
        check_cuda_error(cudaMalloc((void **)&mapping_device, ode_solver->extra_data_size));
        check_cuda_error(cudaMemcpy(mapping_device, mapping, ode_solver->extra_data_size, cudaMemcpyHostToDevice));
    }
    else 
    {
        log_error_and_exit("You need to specify a mask function when using a mixed model!\n");
    }

    solve_gpu<<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps, mapping_device);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    // check_cuda_error(cudaFree(NULL));
    // check_cuda_error(cudaFree(extra_parameters_device));

    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));
    if(mapping_device) check_cuda_error(cudaFree(mapping_device));

    // if(deallocate) free(fibs);
}

__global__ void kernel_set_model_inital_conditions(real *sv, real*IC, int num_volumes)
{
    // Thread ID
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if(threadID < num_volumes) {

        *((real *) ((char *) sv + pitch * 0) + threadID) = 0.f;        // u
        *((real *) ((char *) sv + pitch * 1) + threadID) = 1.f;        // v
        *((real *) ((char *) sv + pitch * 2) + threadID) = 1.f;        // w
        *((real *) ((char *) sv + pitch * 3) + threadID) = 0.f;        // s
    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps, uint32_t *mapping)
{
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;

    // Each thread solves one cell model
    if(threadID < num_cells_to_solve) {
        if(cells_to_solve)
            sv_id = cells_to_solve[threadID];
        else
            sv_id = threadID;

        real rDY[NEQ];

        for (int n = 0; n < num_steps; ++n) {

            if (mapping[sv_id] == 0) {
                RHS_gpu(sv, rDY, stim_currents[threadID], sv_id, dt, 0);
            }
            else if(mapping[sv_id] == 1){
                RHS_gpu(sv, rDY, stim_currents[threadID], sv_id, dt, 1);
            }
            else {
                RHS_gpu(sv, rDY, stim_currents[threadID], sv_id, dt, 2);
            }

            *((real*)((char*)sv) + sv_id) = dt*rDY[0] + *((real*)((char*)sv) + sv_id);
            // *((real*)((char*)sv) + sv_id) = dt*rDY[1] + *((real*)((char*)sv) + sv_id);
            // *((real*)((char*)sv) + sv_id) = dt*rDY[2] + *((real*)((char*)sv) + sv_id);
            // *((real*)((char*)sv) + sv_id) = dt*rDY[3] + *((real*)((char*)sv) + sv_id);

            for(int i = 1; i < NEQ; i++) {
                *((real*)((char*)sv + pitch * i) + sv_id) = dt*rDY[i] + *((real*)((char*)sv + pitch * i) + sv_id);
            }

        }

    }
}


inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_, real dt, int type_cell) {

    const real u   = *((real*)((char*)sv_ + pitch * 0) + threadID_);
    const real v   = *((real*)((char*)sv_ + pitch * 1) + threadID_);
    const real w   = *((real*)((char*)sv_ + pitch * 2) + threadID_);
    const real s   = *((real*)((char*)sv_ + pitch * 3) + threadID_);

    // Constants
    const real u_o = 0.0;
    const real theta_v = 0.3;
    const real theta_w = 0.13;
    const real tau_vplus = 1.4506;
    const real tau_s1 = 2.7342;
    const real k_s = 2.0994;
    const real u_s = 0.9087;

    //real u_o = 0.0;
    real u_u;
    //real theta_v;
    //real theta_w;
    real theta_vminus;
    real theta_o;
    real tau_v1minus;
    real tau_v2minus;
    //real tau_vplus;
    real tau_w1minus;
    real tau_w2minus;
    real k_wminus;
    real u_wminus;
    real tau_wplus;
    real tau_fi;
    real tau_o1;
    real tau_o2;
    real tau_so1;
    real tau_so2;
    real k_so;
    real u_so;
    //real tau_s1;
    real tau_s2;
    //real k_s;
    //real u_s;
    real tau_si;
    real tau_winf;
    real w_infstar;

    if (type_cell == 0) {        // ENDO
      //u_o = 0.0;
      u_u = 1.56;
      //theta_v = 0.3;
      //theta_w = 0.13;
      theta_vminus = 0.2;
      theta_o = 0.006;
      tau_v1minus = 75.0;
      tau_v2minus = 10.0;
      //tau_vplus = 1.4506;
      tau_w1minus = 6.0;
      tau_w2minus = 140.0;
      k_wminus = 200.0;
      u_wminus = 0.016;
      tau_wplus = 280.0;
      tau_fi = 0.1;
      tau_o1 = 470.0;
      tau_o2 = 6.0;
      tau_so1 = 40.0;
      tau_so2 = 1.2;
      k_so = 2.0;
      u_so = 0.65;
      //tau_s1 = 2.7342;
      tau_s2 = 2.0;
      //k_s = 2.0994;
      //u_s = 0.9087;
      tau_si = 2.9013;
      tau_winf = 0.0273;
      w_infstar = 0.78;
    }
    else if (type_cell == 1) {   // MYO
      //u_o = 0.0;
      u_u = 1.61;
      //theta_v = 0.3;
      //theta_w = 0.13;
      theta_vminus = 0.1;
      theta_o = 0.005;
      tau_v1minus = 80.0;
      tau_v2minus = 1.4506;
      //tau_vplus = 1.4506;
      tau_w1minus = 70.0;
      tau_w2minus = 8.0;
      k_wminus = 200.0;
      u_wminus = 0.016;
      tau_wplus = 280.0;
      tau_fi = 0.078;
      tau_o1 = 410.0;
      tau_o2 = 7.0;
      tau_so1 = 91.0;
      tau_so2 = 0.8;
      k_so = 2.1;
      u_so = 0.6;
      //tau_s1 = 2.7342;
      tau_s2 = 4.0;
      //k_s = 2.0994;
      //u_s = 0.9087;
      tau_si = 3.3849;
      tau_winf = 0.01;
      w_infstar = 0.5;
    }
    else {                       // EPI
      //u_o = 0.0;
      u_u = 1.55;
      //theta_v = 0.3;
      //theta_w = 0.13;
      theta_vminus = 0.006;
      theta_o = 0.006;
      tau_v1minus = 60.0;
      tau_v2minus = 1150.0;
      //tau_vplus = 1.4506;
      tau_w1minus = 60.0;
      tau_w2minus = 15.0;
      k_wminus = 65.0;
      u_wminus = 0.03;
      tau_wplus = 200.0;
      tau_fi = 0.11;
      tau_o1 = 400.0;
      tau_o2 = 6.0;
      tau_so1 = 30.0181;
      tau_so2 = 0.9957;
      k_so = 2.0458;
      u_so = 0.65;
      //tau_s1 = 2.7342;
      tau_s2 = 16.0;
      //k_s = 2.0994;
      //u_s = 0.9087;
      tau_si = 1.8875;
      tau_winf = 0.07;
      w_infstar = 0.94;
    }

        // Get du_dt, dv_dt, dw_dt and ds_dt
        //du_dt = - reaction_u(ustep, vstep, wstep, sstep) + I_stim;
        //
        real H = 0.0;
        real h_o = 0.0;
        real h_w = 0.0;
        real h_v_minus = 0.0;
        real tau_vminus = 0.0;
        real J_fi = 0.0;
        real J_so = 0.0;
        real J_si = 0.0;
        real tau_o = 0.0;
        real tau_so = 0.0;
        real tau_s = 0.0;
        // real du_dt = 0.0;
        // real dv_dt = 0.0;
        // real dw_dt = 0.0;
        // real ds_dt = 0.0;
        real v_inf = 0.0;
        real w_inf = 0.0;
        real tau_wminus = 0.0;
        real I_stim = stim_current;

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
        rDY_[0] = - (J_fi + J_so + J_si) + I_stim;
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
        rDY_[1] = (1.0 - H) * (v_inf - v) / tau_vminus - (H * v / tau_vplus);
        //

        //dw_dt = dwdt(ustep, wstep);
        //
        w_inf = (1.0 - h_o) * (1.0 - (u/tau_winf)) + (h_o * w_infstar);
        tau_wminus = tau_w1minus + (((tau_w2minus - tau_w1minus) * (1.0 + tanh(k_wminus*(u - u_wminus)))) * 0.5);
        rDY_[2] = (1.0 - h_w) * (w_inf - w) / tau_wminus - (h_w * w / tau_wplus);
        //

        //ds_dt = dsdt(ustep, sstep);
        //
        tau_s = (1.0 - h_w) * tau_s1 + (h_w * tau_s2);
        rDY_[3] = (((1.0 + tanh(k_s*(u - u_s))) * 0.5) - s) / tau_s;
        //


}
