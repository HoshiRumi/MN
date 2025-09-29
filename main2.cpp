#include <iostream>
#include <cmath>
#include <iomanip>



class NewtonSolver {
private:
    double a, b;
    double eps;
    int max_iter;
    int iterations;

    double f(double x) {
        return x*x*x - 4*x*x + 5*x - 3;
    }


    double df(double x) {
        return 3*x*x - 8*x + 5;
    }


    double d2f(double x) {
        return 6*x - 8;
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
        iterations = 0;

        double x = (f(a) * d2f(a) > 0) ? a : b;
        double x_new;

        std::cout << std::fixed << std::setprecision(6);
        // cout << "\nIterations of Newton's method:\n";
        // cout << "i\t x_i\t\t f(x_i)\t\t |x_{i+1}-x_i|\n";

        while (iterations < max_iter) {
            double fx = f(x);
            double dfx = df(x);

            if (dfx == 0) {
                std::cerr << "Derivative is zero, cannot continue Newton's method." << std::endl;
                return x;
            }

            x_new = x - fx / dfx;
            double diff = fabs(x_new - x);

            // Вывод таблицы (по желанию)
            // cout << iterations+1 << "\t " << x << "\t " << fx << "\t " << diff << "\n";

            iterations++;

            if (diff < eps) break;

            x = x_new;
        }

        return x_new;
    }
};

int main() {
    NewtonSolver solver;
    double a, b, eps;
    char choice;

start:

    std::cout << "Enter the left border a: ";
    std::cin >> a;
    std::cout << "Enter the right border b: ";
    std::cin >> b;
    std::cout << "Enter eps accuracy: ";
    std::cin >> eps;

    solver.setInterval(a, b);
    solver.setEpsilon(eps);
    solver.setMaxIter(1000);

    double root = solver.solve();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nRoot found: x = " << root << std::endl;
    std::cout << "f(x) = " << (root*root*root - 4*root*root + 5*root - 3) << std::endl;
    std::cout << "Iterations: " << solver.getIterations() << std::endl;

    std::cout << "\nWould you like to solve another equation? (y/n): ";
    std::cin >> choice;

    if (choice == 'y' || choice == 'Y') {
        goto start;
    }

    std::cout << "Go read python books" << std::endl;
    return 0;
}
