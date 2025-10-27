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
    std::vector<double> B;                // Вектор правых частей
    int n;                                // Количество уравнений
    int m;                                // Количество неизвестных

    // Функция для вывода расширенной матрицы [A|B]
    void printAugmentedMatrix(const std::vector<std::vector<double>>& augmented, const std::string& title = "") {
        if (!title.empty()) {
            std::cout << "\n" << title << std::endl;
        }
        for (int i = 0; i < augmented.size(); i++) {
            std::cout << "  [ ";
            for (int j = 0; j < augmented[i].size() - 1; j++) {
                std::cout << std::setw(10) << std::fixed << std::setprecision(4) << augmented[i][j] << " ";
            }
            std::cout << "| " << std::setw(10) << std::fixed << std::setprecision(4)
                      << augmented[i][augmented[i].size() - 1] << " ]" << std::endl;
        }
    }

public:
    // Конструктор класса
    LinearEquationsSolver(const std::vector<std::vector<double>>& matrix,
                         const std::vector<double>& vector,
                         int equations,
                         int unknowns) {
        A = matrix;
        B = vector;
        n = equations;
        m = unknowns;
    }

    // Метод Гаусса для решения системы линейных уравнений
    void gaussMethod() {
        // Создаем расширенную матрицу [A|B]
        std::vector<std::vector<double>> augmented(n, std::vector<double>(m + 1));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                augmented[i][j] = A[i][j];
            }
            augmented[i][m] = B[i];
        }

        printAugmentedMatrix(augmented, "Initial augmented matrix:");

        int rank = 0;                    // Текущий ранг матрицы
        std::vector<int> pivotCols;      // Столбцы с опорными элементами
        int stepNumber = 1;

        // Прямой ход метода Гаусса - приведение к ступенчатому виду
        for (int col = 0; col < m && rank < n; col++) {
            // Находим строку с максимальным элементом в текущем столбце (частичный выбор главного элемента)
            int maxRow = rank;
            for (int i = rank + 1; i < n; i++) {
                if (std::abs(augmented[i][col]) > std::abs(augmented[maxRow][col])) {
                    maxRow = i;
                }
            }

            // Если все элементы в столбце близки к нулю, переходим к следующему столбцу
            if (std::abs(augmented[maxRow][col]) < 1e-10) {
                continue;
            }

            // Меняем строки местами, если необходимо
            if (maxRow != rank) {
                std::swap(augmented[rank], augmented[maxRow]);
                std::cout << "\nStep " << stepNumber++ << ": Swap row " << (rank + 1)
                          << " with row " << (maxRow + 1) << std::endl;
                printAugmentedMatrix(augmented, "After row swap:");
            }

            pivotCols.push_back(col);

            // Обнуляем элементы под опорным элементом
            bool hadElimination = false;
            for (int i = rank + 1; i < n; i++) {
                if (std::abs(augmented[i][col]) > 1e-10) {
                    double factor = augmented[i][col] / augmented[rank][col];
                    for (int j = col; j <= m; j++) {
                        augmented[i][j] -= factor * augmented[rank][j];
                    }
                    hadElimination = true;
                }
            }

            if (hadElimination) {
                std::cout << "\nStep " << stepNumber++ << ": Eliminate column " << (col + 1)
                          << " using row " << (rank + 1) << " as pivot" << std::endl;
                printAugmentedMatrix(augmented, "After elimination:");
            }

            rank++;
        }

        std::cout << "\n" << std::string(50, '-') << std::endl;
        printAugmentedMatrix(augmented, "Final row echelon form:");
        std::cout << std::string(50, '-') << std::endl;

        // Проверка на несовместность системы (строки вида 0 = c, где c != 0)
        for (int i = rank; i < n; i++) {
            if (std::abs(augmented[i][m]) > 1e-10) {
                std::cout << "\nSystem is INCONSISTENT (no solution)" << std::endl;
                std::cout << "Reason: Row " << (i + 1) << " gives: 0 = "
                          << std::fixed << std::setprecision(4) << augmented[i][m] << std::endl;
                return;
            }
        }

        // Если ранг меньше числа неизвестных - бесконечно много решений
        if (rank < m) {
            std::cout << "\nSystem has INFINITELY MANY solutions" << std::endl;
            std::cout << "\nGeneral solution:" << std::endl;

            // Определяем свободные переменные
            std::vector<bool> isFreeVar(m, true);
            for (int p : pivotCols) {
                isFreeVar[p] = false;
            }

            // Выражаем базисные переменные через свободные
            for (int i = 0; i < m; i++) {
                if (isFreeVar[i]) {
                    std::cout << "  x" << (i + 1) << " = t" << (i + 1) << " (free parameter)" << std::endl;
                } else {
                    int row = -1;
                    for (int r = 0; r < rank; r++) {
                        if (pivotCols[r] == i) {
                            row = r;
                            break;
                        }
                    }

                    if (row != -1) {
                        std::cout << "  x" << (i + 1) << " = " << std::fixed << std::setprecision(4)
                                  << augmented[row][m] / augmented[row][i];

                        for (int j = 0; j < m; j++) {
                            if (j != i && std::abs(augmented[row][j]) > 1e-10) {
                                double coef = -augmented[row][j] / augmented[row][i];
                                if (coef > 0) {
                                    std::cout << " + " << coef << "*x" << (j + 1);
                                } else {
                                    std::cout << " - " << -coef << "*x" << (j + 1);
                                }
                            }
                        }
                        std::cout << std::endl;
                    }
                }
            }
        } else {
            // Единственное решение - обратный ход метода Гаусса
            std::cout << "\nSystem has UNIQUE solution" << std::endl;
            std::cout << "\nSolution:" << std::endl;

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
                std::cout << "  x" << (i + 1) << " = " << std::fixed << std::setprecision(6) << x[i] << std::endl;
            }
        }
    }

    // Проверка матрицы на диагональное преобладание
    // (|a_ii| > sum(|a_ij|) для всех i, где j != i)
    bool checkDiagonalDominance(const std::vector<std::vector<double>>& matrix) {
        for (int i = 0; i < n; i++) {
            double diagonal = std::abs(matrix[i][i]);
            double sumRow = 0.0;
            for (int j = 0; j < m; j++) {
                if (j != i) {
                    sumRow += std::abs(matrix[i][j]);
                }
            }
            if (diagonal <= sumRow) {
                return false;
            }
        }
        return true;
    }

    // Метод Зейделя (итерационный метод) для решения СЛАУ
    std::pair<std::vector<double>, int> seidelMethod(int maxIterations = 100, double eps = 1e-4) {
        // Метод работает только для квадратных систем
        if (n != m) {
            throw std::runtime_error("Seidel method only works for square systems (n = m)");
        }

        // Проверяем исходную матрицу на диагональное преобладание
        bool isDominant = checkDiagonalDominance(A);

        std::vector<std::vector<double>> A_work = A;
        std::vector<double> B_work = B;

        // Если нет диагонального преобладания, пытаемся переставить строки
        if (!isDominant) {
            std::cout << "\nOriginal matrix does NOT have diagonal dominance." << std::endl;
            std::cout << "Attempting row permutation...\n" << std::endl;

            // Для каждой строки находим максимальный по модулю элемент
            std::vector<int> maxColInRow(n);
            std::vector<double> maxValInRow(n);

            for (int i = 0; i < n; i++) {
                int maxCol = 0;
                double maxVal = std::abs(A[i][0]);
                for (int j = 1; j < m; j++) {
                    if (std::abs(A[i][j]) > maxVal) {
                        maxVal = std::abs(A[i][j]);
                        maxCol = j;
                    }
                }
                maxColInRow[i] = maxCol;
                maxValInRow[i] = maxVal;
            }

            // Пытаемся назначить строки так, чтобы максимальный элемент был на диагонали
            std::vector<int> newRowOrder(n, -1);
            std::vector<bool> usedRows(n, false);

            // Для каждой позиции (столбца диагонали) ищем строку с максимумом в этом столбце
            for (int col = 0; col < m; col++) {
                int bestRow = -1;
                double bestVal = 0;

                for (int row = 0; row < n; row++) {
                    if (usedRows[row]) continue;

                    if (maxColInRow[row] == col) {
                        if (std::abs(A[row][col]) > bestVal) {
                            bestVal = std::abs(A[row][col]);
                            bestRow = row;
                        }
                    }
                }

                // Если не нашли строку с максимумом в этом столбце, берем любую подходящую
                if (bestRow == -1) {
                    for (int row = 0; row < n; row++) {
                        if (usedRows[row]) continue;
                        if (std::abs(A[row][col]) > 1e-10) {
                            if (std::abs(A[row][col]) > bestVal) {
                                bestVal = std::abs(A[row][col]);
                                bestRow = row;
                            }
                        }
                    }
                }

                if (bestRow == -1) {
                    throw std::runtime_error("Cannot apply Seidel method: unable to arrange matrix with diagonal dominance");
                }

                newRowOrder[col] = bestRow;
                usedRows[bestRow] = true;
            }

            // Применяем перестановку строк
            for (int i = 0; i < n; i++) {
                A_work[i] = A[newRowOrder[i]];
                B_work[i] = B[newRowOrder[i]];
            }

            // Проверяем, получилось ли диагональное преобладание после перестановки
            if (!checkDiagonalDominance(A_work)) {
                throw std::runtime_error("Cannot apply Seidel method: unable to achieve diagonal dominance after permutation");
            }

            std::cout << "Row permutation successful!\n" << std::endl;
        } else {
            std::cout << "\nOriginal matrix already has diagonal dominance.\n" << std::endl;
        }

        // Проверка диагональных элементов на ноль
        for (int i = 0; i < m; i++) {
            if (std::abs(A_work[i][i]) < 1e-10) {
                throw std::runtime_error("Diagonal element is zero! Cannot apply Seidel method.");
            }
        }

        // Создаем расширенную матрицу для вывода
        std::vector<std::vector<double>> augmented(n, std::vector<double>(m + 1));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                augmented[i][j] = A_work[i][j];
            }
            augmented[i][m] = B_work[i];
        }

        printAugmentedMatrix(augmented, "Augmented matrix prepared for Seidel method:");

        // Итерационный процесс
        std::vector<double> x(m, 0.0);      // Текущее приближение
        std::vector<double> xNew(m, 0.0);   // Следующее приближение
        int iterations = 0;

        for (int iter = 0; iter < maxIterations; iter++) {
            iterations++;

            // Вычисляем новое приближение
            // x_i^(k+1) = (b_i - sum(a_ij * x_j^(k+1), j<i) - sum(a_ij * x_j^(k), j>i)) / a_ii
            for (int i = 0; i < m; i++) {
                double sum = 0.0;
                for (int j = 0; j < m; j++) {
                    if (j != i) {
                        if (j < i) {
                            sum += A_work[i][j] * xNew[j];  // Используем уже обновленные значения
                        } else {
                            sum += A_work[i][j] * x[j];     // Используем старые значения
                        }
                    }
                }
                xNew[i] = (B_work[i] - sum) / A_work[i][i];
            }

            // Проверяем критерий остановки: max|x_i^(k+1) - x_i^(k)| < eps
            double maxDiff = 0.0;
            for (int i = 0; i < m; i++) {
                maxDiff = std::max(maxDiff, std::abs(xNew[i] - x[i]));
            }

            if (maxDiff < eps) {
                return std::make_pair(xNew, iterations);
            }

            x = xNew;
        }

        return std::make_pair(xNew, iterations);
    }

    // Вывод результатов обоих методов
    void printResults() {
        std::cout << "\n" << std::string(50, '-') << std::endl;
        std::cout << "          GAUSS METHOD (EXACT)" << std::endl;
        std::cout << std::string(50, '-') << std::endl;

        try {
            gaussMethod();
        } catch (const std::exception& e) {
            std::cout << "Error: " << e.what() << std::endl;
        }

        // Метод Зейделя применяется только для квадратных систем
        std::cout << "\n\n" << std::string(50, '-') << std::endl;
        std::cout << "      SEIDEL METHOD (ITERATIVE)" << std::endl;
        std::cout << std::string(50, '-') << std::endl;

        if (n == m) {
            try {
                std::pair<std::vector<double>, int> result = seidelMethod();
                std::vector<double> seidelSolution = result.first;
                int iterations = result.second;

                std::cout << "\nFINAL SOLUTION:" << std::endl;
                for (int i = 0; i < m; i++) {
                    std::cout << "  x" << (i + 1) << " = " << std::fixed
                              << std::setprecision(6) << seidelSolution[i] << std::endl;
                }
                std::cout << "\nTotal iterations: " << iterations << std::endl;
            } catch (const std::exception& e) {
                std::cout << "\nError: " << e.what() << std::endl;
            }
        } else {
            std::cout << "\nSeidel method requires square systems (n = m)" << std::endl;
            std::cout << "Current system: " << n << " equations, " << m << " unknowns" << std::endl;
        }
    }
};

// Основная функция для ввода и решения системы
void solveSystem() {
    int n, m;

    // Ввод размерности системы
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

    // Ввод матрицы коэффициентов A
    std::cout << "\n" << std::string(50, '-') << std::endl;
    std::cout << "Enter coefficient matrix A (" << n << "x" << m << "):" << std::endl;
    std::cout << std::string(50, '-') << std::endl;

    for (int i = 0; i < n; i++) {
        std::cout << "Row " << (i + 1) << " (enter " << m << " coefficients): ";
        for (int j = 0; j < m; j++) {
            std::cin >> matrix[i][j];
        }

        if (std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(10000, '\n');
            std::cout << "Input error! Try again." << std::endl;
            i--;
            continue;
        }
    }

    // Ввод вектора правых частей B
    std::cout << "\n" << std::string(50, '-') << std::endl;
    std::cout << "Enter right-hand side vector B (" << n << " values):" << std::endl;
    std::cout << std::string(50, '-') << std::endl;

    for (int i = 0; i < n; i++) {
        std::cout << "b[" << (i + 1) << "]: ";
        std::cin >> vector[i];

        if (std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(10000, '\n');
            std::cout << "Input error! Try again." << std::endl;
            i--;
            continue;
        }
    }

    // Вывод введенной системы уравнений
    std::cout << "\n" << std::string(70, '-') << std::endl;
    std::cout << "Entered system of equations:" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << "  ";
        for (int j = 0; j < m; j++) {
            if (j > 0 && matrix[i][j] >= 0) std::cout << "+ ";
            std::cout << std::fixed << std::setprecision(2) << matrix[i][j]
                      << "*x" << (j + 1) << " ";
        }
        std::cout << "= " << vector[i] << std::endl;
    }
    std::cout << std::string(70, '-') << std::endl;

    // Создаем объект решателя и запускаем решение
    LinearEquationsSolver solver(matrix, vector, n, m);
    solver.printResults();
}

int main() {
    std::cout << "LINEAR EQUATIONS SOLVER" << std::endl;
    std::cout << "Supports Gauss Method and Seidel Method\n" << std::endl;

    while (true) {
        solveSystem();

        // Меню для продолжения работы или выхода
        std::cout << "\n\n" << std::string(50, '-') << std::endl;
        std::cout << "  What would you like to do next?" << std::endl;
        std::cout << std::string(50, '-') << std::endl;
        std::cout << "  1 - Solve another system" << std::endl;
        std::cout << "  0 - Exit program" << std::endl;
        std::cout << "\nYour choice: ";

        int choice;
        std::cin >> choice;

        if (std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(10000, '\n');
            std::cout << "Invalid input! Please enter 0 or 1." << std::endl;
            continue;
        }

        if (choice == 0) {
            std::cout << "\nThank you for using the program! Goodbye!" << std::endl;
            break;
        } else if (choice != 1) {
            std::cout << "Invalid choice! Please enter 0 or 1." << std::endl;
        }

        std::cout << "\n" << std::string(70, '=') << std::endl;
    }

    return 0;
}
