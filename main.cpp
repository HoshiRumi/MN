#include <iostream>
#include <cmath>
#include <iomanip>

class BisectionSolver {
private:
    double a, b;
    double eps;
    int max_iter;
    int iterations;

    double f(double x) {
        return x*x*x - 4*x*x + 5*x - 3;  // пример: x^3 - 4x^2 + 5x - 3 = 0
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

    double solve() {
        if (f(a) * f(b) >= 0) {
            std::cerr << "Error: the function does not change the sign on ["
                      << a << ", " << b << "]" << std::endl;
            iterations = 0;
            return NAN;
        }

        double c;
        iterations = 0;


        std::cout << "\nIter   a           b           c           f(c)\n";

        while ((b - a) / 2.0 > eps && iterations < max_iter) {
            c = (a + b) / 2.0;


            std::cout << iterations + 1 << "   "
                      << std::fixed << std::setprecision(6)
                      << a << "   "
                      << b << "   "
                      << c << "   "
                      << f(c) << "\n";

            if (f(a) * f(c) < 0) {
                b = c;
            } else {
                a = c;
            }
            iterations++;
        }
        return (a + b) / 2.0;
    }
};

int main() {
    BisectionSolver solver;
    double left, right, eps;
    char choice;

start:

    std::cout << "Enter the left border  a: ";
    std::cin >> left;
    std::cout << "Enter the right border b: ";
    std::cin >> right;
    std::cout << "Enter eps accuracy: ";
    std::cin >> eps;

    solver.setInterval(left, right);
    solver.setEpsilon(eps);
    solver.setMaxIter(1000);

    double root = solver.solve();

    if (!std::isnan(root)) {
        std::cout << "\nRoot found: x = " << std::fixed << std::setprecision(6) << root << std::endl;
        std::cout << "f(x) = " << (root*root*root - 4*root*root + 5*root - 3) << std::endl;
        std::cout << "Iterations: " << solver.getIterations() << std::endl;
    }

    std::cout << "\nWould you like to solve another equation? (y/n): ";
    std::cin >> choice;

    if (choice == 'y' || choice == 'Y') {
        goto start;
    }

    std::cout << "Go read python books" << std::endl;
    return 0;
}
