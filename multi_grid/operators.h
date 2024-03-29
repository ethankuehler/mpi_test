#ifndef C_SOR_3D_OPERATORS_H
#define C_SOR_3D_OPERATORS_H

#include "multi_grid.h"
#include <stdbool.h>
// Takes a N_len can calculates the N_len for a coarsened grid, smallest length must be greater then three.
N_len coarsen(N_len Nlen);

// returns true if grid can be coarsened
bool can_coarsen(N_len Nlen);



//f_in is the fine grid, f_out is the coarse grid, Nlen is that of the fine gird.
//take a fine grind and avenges it to the coarse gird.
void reduce(const double* f_in, double* f_out, N_len Nlen);

/*
 * calculates the residual of u.
 * f is the right hand side of the dif eq.
 * u is the approximate solution of the dif eq.
 * r is the output.
 * Nlen is the dimension of the gird.
 * dxs is the step size squared.
 */
void residual(const double* f, const double* u, double* r, N_len Nlen, double dxs);

/*
 * This function calculates the residual of u and then reducing the residual to an m-1 grid
 * f is the right hand side of the dif eq
 * u is the approximate solution
 * f_out is the output
 * Nlen is the dimension of the fine gird.
 * dxs is the step size squared
 */
void restriction(const double* f, const double* u, double* f_out, N_len Nlen, double dxs);

/*
 * interpolates the function f into f_out where m is that of the coarse grid.
 * Nlen is that of the finer grid.
 * dx is the step size.
 */
void interpolate(const double* f, double* f_out, N_len Nlen, double dx);

#endif //C_SOR_3D_OPERATORS_H
