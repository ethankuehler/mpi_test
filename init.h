#ifndef MPI_TEST_INIT_H
#define MPI_TEST_INIT_H

#include "multi_grid/multi_grid.h"
#include "tSolver/tSolver.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct _com {
    bool read;
    int rank;
}Com;

typedef struct grid_ {
    N_len Nlen;
    N_len position;
    double* data;
    int rank;
}Grid;

void init_com(Com* com, int rank);

void init_source(double* f, N_len Nlen, N_len position, N_len middle, double dx, double dens);
void init_guess(double* u, const double* f, N_len Nlen, N_len middle, double dx, double shift, double dens);

/*
 * mu or master u is the larger gird that is made up of all other girds while u is a sub gird
 */
void combine_grid(double* mu, const double* u, N_len mu_Nlen, N_len u_Nlen, N_len pos);

void send_grid(Grid* grid, int receiver);

void recv_grid(Grid* grid, int sender);

void save_gird(const char* str, double* vec, N_len Nlen);


#endif //MPI_TEST_INIT_H
