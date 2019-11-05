#include <stdio.h>
#include "multi_grid/multi_grid.h"
#include "tSolver/tSolver.h"
#include "init.h"
#include <mpi.h>
#include <stdlib.h>
void solve_top(const double* f, double* u, N_len Nlen, double dx) {
    tSolve(f, u, Nlen, 5, 1.9, dx);
}

void solve_coarse(const double* f, double* u, N_len Nlen, double dx) {
    tSolve(f, u, Nlen, 5, 1, dx);
}

void solve_base(const double* f, double* u, N_len Nlen, double dx) {
    tSolve(f, u, Nlen, 1, 1, dx);
}

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);


    int N = 128;

    double dx = 0.1;

    int L = 2*N;
    N_len middle = (N_len){L / 2, L / 2, L / 2};

    Grid grid;
    grid.Nlen = (N_len){2*N, N, N};
    grid.data = calloc(sizeof(double), length(grid.Nlen));
    grid.rank = rank;

    if(rank == 0) {
        grid.position = (N_len){0, 0, 0};
    } else if(rank == 1) {
        grid.position = (N_len){0, 0, N};
    } else if(rank == 2) {
        grid.position = (N_len){0, N, 0};
    } else if(rank == 3) {
        grid.position = (N_len){0, N, N};
    }

    init_source(grid.data, grid.Nlen, grid.position, middle, dx, 1);


    if (rank != 0) {
        send_grid(&grid, 0);
    } else {

        Grid* grids = calloc(sizeof(Grid), world_size);
        for(int i = 1; i < world_size; i++) {
            recv_grid(&grids[i], i);
        }

        grids[0] = grid;
        Grid main_grid;
        main_grid.Nlen = (N_len){2*N + 1, 2*N + 1, 2*N + 1};
        main_grid.position = (N_len){0, 0, 0};
        main_grid.data = calloc(sizeof(double), length(main_grid.Nlen));
        main_grid.rank = 0;
        int i = 0;
#pragma omp parallel for num_threads(4)
        for( i = 0; i < world_size; i++) {
            combine_grid(main_grid.data, grids[i].data, main_grid.Nlen, grids[i].Nlen, grids[i].position);
        }

        save_gird("data4", main_grid.data, main_grid.Nlen);
        Grid solve = main_grid;
        solve.data = calloc(sizeof(double), length(solve.Nlen));

        init_guess(solve.data, main_grid.data, solve.Nlen, middle, dx, 0.1, 1);
        funcs_args arg = (funcs_args){solve_top, solve_coarse, solve_base};
        multi(main_grid.data, solve.data, solve.Nlen, dx, arg, true);

        save_gird("data5", solve.data, solve.Nlen);

    }



    char* str = calloc(sizeof(char), 10);
    sprintf(str, "data%d", rank);
    save_gird(str, grid.data, grid.Nlen);

    // Finalize the MPI environment.
    MPI_Finalize();
}