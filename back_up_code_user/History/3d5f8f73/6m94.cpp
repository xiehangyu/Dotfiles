#include <iostream>
#include <vector>
#include <cmath>

void doolittleLU(std::vector<std::vector<double>>& A) {
    int n = A.size();

    for (int i = 0; i < n; ++i) {
        // 计算U
        for (int j = i; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < i; ++k) {
                sum += A[i][k] * A[k][j];
            }
            A[i][j] = A[i][j] - sum;
        }

        // 计算L
        for (int j = i + 1; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < i; ++k) {
                sum += A[j][k] * A[k][i];
            }
            A[j][i] = (A[j][i] - sum) / A[i][i];
        }
    }
}

int main() {
    std::vector<std::vector<double>> A = {
        {4, 3, 2, 1},
        {5, 4, 3, 2},
        {6, 5, 4, 3},
        {7, 6, 5, 4}
    };

    doolittleLU(A);

    // 打印分解后的矩阵
    for (const auto &row : A) {
        for (const auto &elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

void swapRows(std::vector<std::vector<double>>& A, int row1, int row2) {
    std::swap(A[row1], A[row2]);
}

void doolittleLUWithPartialPivoting(std::vector<std::vector<double>>& A) {
    int n = A.size();

    for (int i = 0; i < n; ++i) {
        // 部分主元选取
        int maxIndex = i;
        double maxValue = std::abs(A[i][i]);
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > maxValue) {
                maxIndex = k;
                maxValue = std::abs(A[k][i]);
            }
        }

        // 交换行
        if (maxIndex != i) {
            swapRows(A, i, maxIndex);
        }

        // 计算U
        for (int j = i; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < i; ++k) {
                sum += A[i][k] * A[k][j];
            }
            A[i][j] = A[i][j] - sum;
        }

        // 计算L
        for (int j = i + 1; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < i; ++k) {
                sum += A[j][k] * A[k][i];
            }
            A[j][i] = (A[j][i] - sum) / A[i][i];
        }
    }
}

int main() {
    std::vector<std::vector<double>> A = {
        {1, 3, 2, 1},
        {5, 4, 3, 2},
        {6, 5, 4, 3},
        {7, 6, 5, 4}
    };

    doolittleLUWithPartialPivoting(A);

    // 打印分解后的矩阵
    for (const auto &row : A) {
        for (const auto &elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
