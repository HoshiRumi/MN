#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>

class LinearEquationsSolver {
private:
    std::vector<std::vector<double>> A;  // Матрица коэффициентов
    std::vector<double> B;               // Вектор правых частей
    int n;                               // Количество уравнений
    int m;                               // Количество неизвестных

public:
    // Конструктор
    LinearEquationsSolver(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector, int equations, int unknowns) {
        A = matrix;
        B = vector;
        n = equations;
        m = unknowns;
    }

    // Метод Гаусса с определением ранга
    void gaussMethod() {
        // Создаем расширенную матрицу
        std::vector<std::vector<double>> augmented(n, std::vector<double>(m + 1));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                augmented[i][j] = A[i][j];
            }
            augmented[i][m] = B[i];
        }

        // Прямой ход
        int rank = 0;
        std::vector<int> pivotCols;

        for (int col = 0; col < m && rank < n; col++) {
            // Поиск главного элемента
            int maxRow = rank;
            for (int i = rank + 1; i < n; i++) {
                if (std::abs(augmented[i][col]) > std::abs(augmented[maxRow][col])) {
                    maxRow = i;
                }
            }

            if (std::abs(augmented[maxRow][col]) < 1e-10) {
                continue; // Пропускаем этот столбец
            }

            // Перестановка строк
            std::swap(augmented[rank], augmented[maxRow]);
            pivotCols.push_back(col);

            // Обнуление элементов
            for (int i = rank + 1; i < n; i++) {
                double factor = augmented[i][col] / augmented[rank][col];
                for (int j = col; j <= m; j++) {
                    augmented[i][j] -= factor * augmented[rank][j];
                }
            }
            rank++;
        }

        // Проверка совместности
        for (int i = rank; i < n; i++) {
            if (std::abs(augmented[i][m]) > 1e-10) {
                std::cout << "System is inconsistent (no solution)." << std::endl;
                return;
            }
        }

        std::cout << "Matrix rank: " << rank << std::endl;

        if (rank < m) {
            std::cout << "System has infinitely many solutions (underdetermined)." << std::endl;
            std::cout << "Free variables: " << (m - rank) << std::endl;
            std::cout << "\nParticular solution (setting free variables to 0):" << std::endl;
        } else {
            std::cout << "System has unique solution." << std::endl;
            std::cout << "\nSolution:" << std::endl;
        }

        // Обратный ход
        std::vector<double> x(m, 0.0);
        for (int i = rank - 1; i >= 0; i--) {
            int col = pivotCols[i];
            x[col] = augmented[i][m];
            for (int j = col + 1; j < m; j++) {
                x[col] -= augmented[i][j] * x[j];
            }
            x[col] /= augmented[i][col];
        }

        for (int i = 0; i < m; i++) {
            std::cout << "  x" << (i + 1) << " = " << std::fixed << std::setprecision(6) << x[i];
            bool isFree = true;
            for (int p : pivotCols) {
                if (p == i) {
                    isFree = false;
                    break;
                }
            }
            if (isFree && rank < m) {
                std::cout << " (free variable, set to 0)";
            }
            std::cout << std::endl;
        }
    }

    // Метод Зейделя (только для квадратных систем с диагональным преобладанием)!!!
    std::pair<std::vector<double>, int> seidelMethod(int maxIterations = 1000, double tolerance = 1e-6) {
        if (n != m) {
            throw std::runtime_error("Seidel method only works for square systems (n = m)");
        }

        std::vector<double> x(m, 0.0);
        std::vector<double> xNew(m, 0.0);
        int iterations = 0;

        // Проверка условия диагонального преобладания
        for (int i = 0; i < m; i++) {
            double diagonal = std::abs(A[i][i]);
            double sumRow = 0.0;
            for (int j = 0; j < m; j++) {
                if (j != i) {
                    sumRow += std::abs(A[i][j]);
                }
            }
            if (diagonal <= sumRow) {
                std::cout << "Warning: diagonal dominance condition is not satisfied." << std::endl;
                std::cout << "Seidel method may not converge!" << std::endl;
                break;
            }
        }

        // Итерационный процесс
        for (int iter = 0; iter < maxIterations; iter++) {
            iterations++;

            for (int i = 0; i < m; i++) {
                double sum1 = 0.0;
                for (int j = 0; j < i; j++) {
                    sum1 += A[i][j] * xNew[j];
                }

                double sum2 = 0.0;
                for (int j = i + 1; j < m; j++) {
                    sum2 += A[i][j] * x[j];
                }

                xNew[i] = (B[i] - sum1 - sum2) / A[i][i];
            }

            // Проверка критерия остановки
            double maxDiff = 0.0;
            for (int i = 0; i < m; i++) {
                maxDiff = std::max(maxDiff, std::abs(xNew[i] - x[i]));
            }

            if (maxDiff < tolerance) {
                return std::make_pair(xNew, iterations);
            }

            x = xNew;
        }

        std::cout << "Warning: maximum number of iterations reached (" << maxIterations << ")" << std::endl;
        return std::make_pair(xNew, iterations);
    }

    // Вывод результатов
    void printResults() {
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "LINEAR SYSTEM SOLUTION RESULTS" << std::endl;
        std::cout << std::string(60, '=') << std::endl;

        // Метод Гаусса
        std::cout << "\n--- GAUSS METHOD ---" << std::endl;
        try {
            gaussMethod();
        } catch (const std::exception& e) {
            std::cout << "Error: " << e.what() << std::endl;
        }

        // Метод Зейделя (только для квадратных систем)
        if (n == m) {
            std::cout << "\n--- SEIDEL METHOD ---" << std::endl;
            try {
                std::pair<std::vector<double>, int> result = seidelMethod();
                std::vector<double> seidelSolution = result.first;
                int iterations = result.second;

                std::cout << "Solution:" << std::endl;
                for (int i = 0; i < m; i++) {
                    std::cout << "  x" << (i + 1) << " = " << std::fixed << std::setprecision(6) << seidelSolution[i] << std::endl;
                }
                std::cout << "\nNumber of iterations: " << iterations << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
            }
        } else {
            std::cout << "\n--- SEIDEL METHOD ---" << std::endl;
            std::cout << "Seidel method is only applicable to square systems (n = m)." << std::endl;
            std::cout << "Current system: " << n << " equations, " << m << " unknowns." << std::endl;
        }

        std::cout << "\n" << std::string(60, '=') << std::endl;
    }
};

int main() {
    int n, m;

    // Ввод размерности
    while (true) {
        std::cout << "\nEnter number of equations (n): ";
        std::cin >> n;

        if (std::cin.fail() || n <= 0) {
            std::cin.clear();
            std::cin.ignore(10000, '\n');
            std::cout << "Error! Enter a positive integer." << std::endl;
            continue;
        }
        break;
    }

    while (true) {
        std::cout << "Enter number of unknowns (m): ";
        std::cin >> m;

        if (std::cin.fail() || m <= 0) {
            std::cin.clear();
            std::cin.ignore(10000, '\n');
            std::cout << "Error! Enter a positive integer." << std::endl;
            continue;
        }
        break;
    }

    std::vector<std::vector<double>> matrix(n, std::vector<double>(m));
    std::vector<double> vector(n);

    // Ввод коэффициентов
    std::cout << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << "Equation " << (i + 1) << ": ";

        for (int j = 0; j < m; j++) {
            std::cin >> matrix[i][j];
        }
        std::cin >> vector[i];

        if (std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(10000, '\n');
            std::cout << "Input error! Try again." << std::endl;
            i--;
            continue;
        }
    }

    // Вывод введенной системы
    std::cout << "\nEntered system of equations:" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << "  ";
        for (int j = 0; j < m; j++) {
            if (j > 0 && matrix[i][j] >= 0) std::cout << "+ ";
            std::cout << std::fixed << std::setprecision(2) << matrix[i][j] << "*x" << (j + 1) << " ";
        }
        std::cout << "= " << vector[i] << std::endl;
    }

    // Решение системы
    LinearEquationsSolver solver(matrix, vector, n, m);
    solver.printResults();

    return 0;
}
