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


    return 0;
}
