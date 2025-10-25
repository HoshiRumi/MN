#include <iostream>
#include <cmath>
#include <iomanip>

class EquationSolver {
private:
    double a, b;
    double eps;
    int max_iter;
    int iterations;

    // Функция и ее производные
    double f(double x) {
        return 2*x*x*x - 0.6*x*x + 0.6*x - 0.2;
    }

    double df(double x) {
        return 6*x*x - 1.2*x + 0.6;
    }

    double d2f(double x) {
        return 12*x - 1.2;
    }

    // Проверка условия смены знака
    bool checkSignChange() {
        return f(a) * f(b) < 0;
    }

public:
    void setInterval(double left, double right) {
        a = left;
        b = right;
    }

    void setEpsilon(double epsilon) {
        eps = epsilon;
    }

    void setMaxIter(int iter) {
        max_iter = iter;
    }

    int getIterations() const {
        return iterations;
    }

    // Публичный метод для получения значения функции
    double getFunctionValue(double x) {
        return f(x);
    }

    // Метод бисекции
    double solveBisection() {
        if (!checkSignChange()) {
            iterations = 0;
            return NAN;
        }

        double left = a, right = b;
        double c;
        iterations = 0;

        while ((right - left) / 2.0 > eps && iterations < max_iter) {
            c = (left + right) / 2.0;

            if (f(left) * f(c) < 0) {
                right = c;
            } else {
                left = c;
            }
            iterations++;
        }
        return (left + right) / 2.0;
    }

    // Метод Ньютона
    double solveNewton() {
        if (!checkSignChange()) {
            iterations = 0;
            return NAN;
        }

        iterations = 0;
        double x = (f(a) * d2f(a) > 0) ? a : b;
        double x_new;

        while (iterations < max_iter) {
            double fx = f(x);
            double dfx = df(x);

            if (dfx == 0) {
                return x;
            }

            x_new = x - fx / dfx;
            double diff = fabs(x_new - x);

            iterations++;

            if (diff < eps) break;

            x = x_new;
        }

        return x_new;
    }

    // Метод секущих
    double solveSecant() {
        if (!checkSignChange()) {
            iterations = 0;
            return NAN;
        }

        iterations = 0;

        double x0, x1;
        if (f(a) * d2f(a) > 0) {
            x0 = a;
            x1 = a + eps;
        } else {
            x0 = b;
            x1 = b - eps;
        }

        double x_prev = x0;
        double x_curr = x1;
        double x_next = x_curr;

        while (iterations < max_iter) {
            double f_prev = f(x_prev);
            double f_curr = f(x_curr);

            if (f_curr - f_prev == 0) {
                return x_curr;
            }

            x_next = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev);
            double diff = fabs(x_next - x_curr);

            iterations++;

            if (diff < eps) break;

            x_prev = x_curr;
            x_curr = x_next;
        }

        return x_next;
    }
};

int main() {
    EquationSolver solver;
    double left, right, eps;
    char choice;

start:
    std::cout << "Enter the left border a: ";
    std::cin >> left;
    std::cout << "Enter the right border b: ";
    std::cin >> right;
    std::cout << "Enter eps accuracy: ";
    std::cin >> eps;

    solver.setInterval(left, right);
    solver.setEpsilon(eps);
    solver.setMaxIter(1000);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "COMPARISON OF ROOT-FINDING METHODS" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    // Метод бисекции
    std::cout << "\n1. BISECTION METHOD:" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    double root_bisection = solver.solveBisection();
    if (!std::isnan(root_bisection)) {
        std::cout << "Root found: x = " << root_bisection << std::endl;
        std::cout << "f(x) = " << solver.getFunctionValue(root_bisection) << std::endl;
        std::cout << "Iterations: " << solver.getIterations() << std::endl;
    } else {
        std::cout << "Method failed!" << std::endl;
    }

    // Метод Ньютона
    std::cout << "\n2. NEWTON'S METHOD:" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    double root_newton = solver.solveNewton();
    if (!std::isnan(root_newton)) {
        std::cout << "Root found: x = " << root_newton << std::endl;
        std::cout << "f(x) = " << solver.getFunctionValue(root_newton) << std::endl;
        std::cout << "Iterations: " << solver.getIterations() << std::endl;
    } else {
        std::cout << "Method failed!" << std::endl;
    }

    // Метод секущих
    std::cout << "\n3. SECANT METHOD:" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    double root_secant = solver.solveSecant();
    if (!std::isnan(root_secant)) {
        std::cout << "Root found: x = " << root_secant << std::endl;
        std::cout << "f(x) = " << solver.getFunctionValue(root_secant) << std::endl;
        std::cout << "Iterations: " << solver.getIterations() << std::endl;
    } else {
        std::cout << "Method failed!" << std::endl;
    }

    std::cout << "\n" << std::string(60, '=') << std::endl;

    std::cout << "\nWould you like to solve another equation? (y/n): ";
    std::cin >> choice;

    if (choice == 'y' || choice == 'Y') {
        goto start;
    }

    std::cout << "Go read python books" << std::endl;
    return 0;
}
