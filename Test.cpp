#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <unistd.h>

// Подключаем функции из основного кода
double recur(vector<vector<double>>& C, double dT, double h, int n, int j) {
    return C[n][j] - n * dT * dT * (C[n][j] - C[n][j - 1]) / h + 4 * (j - 1) * h * dT;
}

// Тест для проверки начального условия
TEST(InitialConditionTest, CheckInitialCondition) {
    double dT = 0.1;
    double h = 0.1;
    double x_0 = 0;
    double x_n = 1;
    double t_0 = 0;
    double t_n = 0;

    int x_nums = int((x_n - x_0) / h);
    int N = int((t_n - t_0) / dT);

    vector<vector<double>> C(N + 1, vector<double>(x_nums + 1));

    // Проверяем, что начальное условие c(0, x) = sin(x) + 1 правильно установлено
    for (int i = 0; i <= x_nums; i++) {
        ASSERT_DOUBLE_EQ(C[0][i], sin(x_0 + h * i) + 1.0);
    }
}

// Тест для проверки граничного условия
TEST(BorderConditionTest, CheckBorderCondition) {
    double dT = 0.1;
    double h = 0.1;
    double x_0 = 0;
    double x_n = 0;
    double t_0 = 0;
    double t_n = 1;

    int x_nums = int((x_n - x_0) / h);
    int N = int((t_n - t_0) / dT);

    vector<vector<double>> C(N + 1, vector<double>(x_nums + 1));

    // Проверяем, что граничное условие c(t, 0) = 1 правильно установлено
    for (int i = 0; i <= N; i++) {
        ASSERT_DOUBLE_EQ(C[i][0], 1.0);
    }
}

// Тест для проверки размерности матрицы
TEST(MatrixSizeTest, CheckMatrixSize) {
    double dT = 0.1;
    double h = 0.1;
    double x_0 = 0;
    double x_n = 1;
    double t_0 = 0;
    double t_n = 1;

    int x_nums = int((x_n - x_0) / h);
    int N = int((t_n - t_0) / dT);

    vector<vector<double>> C(N + 1, vector<double>(x_nums + 1));


    // Проверяем, что матрица имеет правильные размеры
    ASSERT_EQ(C.size(), N + 1);
    ASSERT_EQ(C[0].size(), x_nums + 1);
}

// Тест для проверки рекуррентного соотношения
TEST(RecurrenceRelationTest, CheckRecurrenceRelation) {
    double dT = 0.1;
    double h = 0.1;
    double x_0 = 0;
    double x_n = 1;
    double t_0 = 0;
    double t_n = 1;

    int x_nums = int((x_n - x_0) / h);
    int N = int((t_n - t_0) / dT);

    vector<vector<double>> C(N + 1, vector<double>(x_nums + 1));

    // Инициализируем начальные условия
    for (int x_index = 0; x_index <= x_nums; x_index++) {
        C[0][x_index] = sin(x_0 + h * x_index) + 1.0;
    }

    // Запускаем рекуррентное соотношение
    for (int n = 0; n < N; n++) {
        for (int x_index = 1; x_index <= x_nums; x_index++) {
            double result = recur(C, dT, h, n, x_index);
            ASSERT_DOUBLE_EQ(C[n + 1][x_index], result);
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
