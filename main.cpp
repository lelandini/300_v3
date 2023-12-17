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

    // Проверка, является ли size делителем x_nums
    if (x_nums % size != 0) {
        if (rank == 0){
            cerr << "Ошибка! Количество процессов не является подходящим" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    int local_size = x_nums / size;
    vector<vector<double>> C(N + 1, vector<double>(x_nums + 1));

    // Установка начальных условий
    for (int x_index = 0; x_index <= x_nums; x_index++) {
        C[0][x_index] = sin(x_0 + h * x_index) + 1.0;
        for (int N_index = 1; N_index <= N; N_index++){
            C[N_index][0] = 1.0;
        }
    }

    // Решение уравнения
    for (int n = 0; n < N; n++) {

        MPI_Bcast(&C[n][0], C[0].size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (int x_index = rank * local_size + 1; x_index <= (rank + 1) * local_size; x_index++) {
            C[n + 1][x_index] = recur(C, dT, h, n, x_index);
        }

        // Сбор результатов на нулевом процессе
        MPI_Gather(&C[n + 1][rank * local_size + 1], local_size, MPI_DOUBLE, &C[n + 1][1], local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    // Вывод результатов на 0 процессе
    if (rank == 0) {
        cout << "dT = " << dT << "; h = " << h << '\n';
        cout << "| t\\x|";
        for (int i = 0; i < x_nums + 1; i++) {
            cout << setw(8) << fixed << setprecision(3) << x_0 + h * i << "|";
        }
        cout << endl;

        for (int i = 0; i < N + 1; i++) {
            cout << setw(5) << fixed << setprecision(3) << t_0 + dT * i << "|";
            for (int j = 0; j < x_nums + 1; j++) {
                cout << setw(8) << fixed << setprecision(3) << C[i][j] << "|";
            }
            cout << endl;
        }
    }

    return 0;
}
