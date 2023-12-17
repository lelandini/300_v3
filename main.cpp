#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <mpi.h>
#include <unistd.h>

using namespace std;

double recur(vector<vector<double>>& C, double dT, double h, int n, int j) {
    return C[n][j] - n * dT * dT * (C[n][j] - C[n][j - 1]) / h + 4 * (j - 1) * h * dT;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double dT = 0.1;
    double h = 0.1;
    double x_0 = 0;
    double x_n = 1;
    double t_0 = 0;
    double t_n = 1;

    int x_nums = int((x_n - x_0) / h);
    int N = int((t_n - t_0) / dT);

    int local_size = x_nums / size;
    vector<vector<double>> C(N + 1, vector<double>(x_nums + 1));

    // Установка начальных условий
    for (int x_index = 0; x_index <= x_nums; x_index++) {
        C[0][x_index] = sin(x_0 + h * x_index) + 1.0;
        for (int N_index = 1; N_index <= N; N_index++){
            C[N_index][0] = 1.0;
        }
    }

    return 0;
}
