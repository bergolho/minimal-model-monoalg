// #ifndef FUNCTIONS_INC_C
// #define FUNCTIONS_INC_C

// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <math.h>
// #include <stdbool.h>
// #include <omp.h>

// /*-----------------------------------------------------
// Auxiliary functions
// -----------------------------------------------------*/
// // Standard Heaviside function
// double H(double x)
// {
//     if (x > 0.0)
//     {
//         return 1.0;
//     }
//     else
//     {
//         return 0.0;
//     }
// }


// /*-----------------------------------------------------
// Functions of voltage variable u
// -----------------------------------------------------*/
// double tau_vminus(double u)
// {
//     double h = H(u - theta_vminus);
//     return (1.0 - h) * tau_v1minus + (h * tau_v2minus);
// }

// double tau_wminus(double u)
// {
//     return tau_w1minus + (((tau_w2minus - tau_w1minus) * (1.0 + tanh(k_wminus*(u - u_wminus)))) * 0.5);
// }

// double tau_so(double u)
// {
//     return tau_so1 + (((tau_so2 - tau_so1) * (1.0 + tanh(k_so*(u - u_so)))) * 0.5);
// }

// double tau_s(double u)
// {
//     double h = H(u - theta_w);
//     return (1.0 - h) * tau_s1 + (h * tau_s2);
// }

// double tau_o(double u)
// {
//     double h = H(u - theta_o);
//     return (1.0 - h) * tau_o1 + (h * tau_o2);
// }

// // Infinity values
// double v_inf_function(double u)
// {
//     if (u < theta_vminus)
//     {
//         return 1.0;
//     }
//     else
//     {
//         return 0.0;
//     }
// }

// double w_inf_function(double u)
// {
//     double h = H(u - theta_o);
//     return (1.0 - h) * (1.0 - (u/tau_winf)) + (h * w_infstar);
// }


// /*-----------------------------------------------------
// Currents functions
// -----------------------------------------------------*/
// double J_fi(double u, double v)
// {
//     return -v * H(u-theta_v) * (u-theta_v) * (u_u-u) / tau_fi;
// }

// double J_so(double u)
// {
//     double h = H(u-theta_w);
//     return ((u-u_o) * (1.0 - h) / tau_o(u)) + (h / tau_so(u));
// }

// double J_si(double u, double w, double s)
// {
//     return - H(u-theta_w) * w * s / tau_si;
// }


// /*-----------------------------------------------------
// Differential equations for each variable
// -----------------------------------------------------*/
// double reaction_u(double u, double v, double w, double s)
// {
//     return (J_fi(u, v) + J_so(u) + J_si(u, w, s));
// }

// double dvdt(double u, double v)
// {
//     double h = H(u - theta_v);
//     return (1.0 - h) * (v_inf_function(u) - v) / tau_vminus(u) - (h * v / tau_vplus);
// }

// double dwdt(double u, double w)
// {
//     double h = H(u - theta_w);
//     return (1.0 - h) * (w_inf_function(u) - w) / tau_wminus(u) - (h * w / tau_wplus);
// }

// double dsdt(double u, double s)
// {
//     return (((1.0 + tanh(k_s*(u - u_s))) * 0.5) - s) / tau_s(u);
// }


// #endif // FUNCTIONS_INC_C