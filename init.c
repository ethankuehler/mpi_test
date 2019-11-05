#include "multi_grid/multi_grid.h"
#include "init.h"
#include <math.h>
#include <mpi.h>

void init_com(Com* com, int rank) {
    com->rank = rank;
    com->read = false;
}

void init_source(double* f, N_len Nlen, N_len position, N_len middle, double dx, double dens) {
    for(int i = 0; i < Nlen.i; i++) {
        for(int j = 0; j < Nlen.j; j++) {
            for(int k = 0; k < Nlen.k; k++) {
                double x = i*dx + position.i*dx - middle.i*dx;
                double y = j*dx + position.j*dx - middle.j*dx;
                double z = k*dx + position.k*dx - middle.k*dx;
                double R = sqrt(x*x + y*y + z*z);
                int n = loc(i, j, k, Nlen);
                if (R < 1) {
                    f[n] = M_PI*4*dens;
                } else {
                    f[n] = 0;
                }
            }
        }
    }
}

void init_guess(double* u, const double* f, N_len Nlen, N_len middle, double dx, double shift, double dens) {

    //integrate mass
    double m = 0;
    for (int i = 0; i < Nlen.i; i++) {
        for (int j = 0; j < Nlen.j; j++) {
            for (int k = 0; k < Nlen.k; k++) {
                int n = loc(i,j,k, Nlen);
                if (f[n] != 0.0) {
                    m += dens*dx*dx*dx;
                }
            }
        }
    }

    for(int i = 0; i < Nlen.i; i++) {
        for (int j = 0; j < Nlen.j; j++) {
            for (int k = 0; k < Nlen.k; k++) {
                double x = i*dx - middle.i*dx;
                double y = j*dx - middle.j*dx;
                double z = k*dx - middle.k*dx + shift;
                double R = 1;
                double rsqrd = x*x + y*y + z*z;
                int n = loc(i, j, k, Nlen);
                if (rsqrd < R*R) {
                    u[n] = -(m/(2*R*R*R))*(3*R*R - rsqrd);
                } else {
                    u[n] = -m/sqrt(rsqrd);
                }
            }
        }
    }
}

void combine_grid(double* mu, const double* u, N_len mu_Nlen, N_len u_Nlen, N_len pos) {
    for(int i = 0; i < u_Nlen.i; i++) {
        for(int j = 0; j < u_Nlen.j; j++) {
            for(int k = 0; k < u_Nlen.k; k++) {
                //position on the master grid
                int x = i + pos.i;
                int y = j + pos.j;
                int z = k + pos.k;
                int mu_n = loc(x,y,z, mu_Nlen);
                int u_n = loc(i,j,k, u_Nlen);
                mu[mu_n] = u[u_n];
            }
        }
    }
}

void send_grid(Grid* grid, int receiver) {
    MPI_Send(&grid->Nlen,3, MPI_INT, receiver,1, MPI_COMM_WORLD);
    MPI_Send(&grid->position,3, MPI_INT, receiver,1, MPI_COMM_WORLD);
    printf("send %i\n", length(grid->Nlen));
    MPI_Send(grid->data,length(grid->Nlen), MPI_DOUBLE, receiver,1, MPI_COMM_WORLD);
    MPI_Send(&grid->rank,1, MPI_INT, receiver,1, MPI_COMM_WORLD);
}

void recv_grid(Grid* grid, int sender) {
    MPI_Recv(&grid->Nlen, 3, MPI_INT, sender, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&grid->position, 3, MPI_INT, sender, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    grid->data = calloc(sizeof(double), length(grid->Nlen));
    printf("recv %i\n", length(grid->Nlen));
    MPI_Recv(grid->data, length(grid->Nlen), MPI_DOUBLE, sender, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&grid->rank, 1, MPI_INT, sender, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void save_gird(const char* str, double* vec, N_len Nlen) {
    int l = length(Nlen);
    char* txt = ".txt";
    char* fmt = ".fmt";
    char* loc1 = calloc(sizeof(char), strlen(str) + 5);
    char* loc2 = calloc(sizeof(char), strlen(str) + 5);
    strcpy(loc1, str);
    strcpy(loc2, str);
    strcat(loc1, txt);
    strcat(loc2, fmt);


    FILE* file = fopen(loc1, "w");

    for (int i = 0; i < l; i++) {
        fprintf(file, "%f ", vec[i]);
    }

    fprintf(file, "\n");
    fclose(file);

    file = fopen(loc2, "w");
    fprintf(file, "%d ", Nlen.i);
    fprintf(file, "%d ", Nlen.j);
    fprintf(file, "%d\n", Nlen.k);
    fclose(file);
}
