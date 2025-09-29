#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

class SecantSolver {
private:
    double a, b;
    double eps;
    int max_iter;
    int iterations;

    double f(double x) {
        return x*x*x - 4*x*x + 5*x - 3;
    }

    double d2f(double x) {
        return 6*x - 8; 
    }

public:
    void setInitialPoints(double left, double right) {
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


        double x0, x1;
        if (f(a) * d2f(a) > 0) {
            x0 = a;
            x1 = b;
        } else {
            x0 = b;
            x1 = a;
        }

        double x_prev = x0;
        double x_curr = x1;
        double x_next;

        while (iterations < max_iter) {
            double f_prev = f(x_prev);
            double f_curr = f(x_curr);

            if (f_curr - f_prev == 0) {
                cerr << "Division by zero in Secant method" << endl;
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
    SecantSolver solver;
    double a, b, eps;
    char choice;

start:

    cout << "Enter the left border a: ";
    cin >> a;
    cout << "Enter the right border b: ";
    cin >> b;
    cout << "Enter eps accuracy: ";
    cin >> eps;

    solver.setInitialPoints(a, b);
    solver.setEpsilon(eps);
    solver.setMaxIter(1000);

    double root = solver.solve();

    cout << fixed << setprecision(6);
    cout << "\nRoot found: x = " << root << endl;
    cout << "f(x) = " << (root*root*root - 4*root*root + 5*root - 3) << endl;
    cout << "Iterations: " << solver.getIterations() << endl;

    cout << "\nWould you like to solve another equation? (y/n): ";
    cin >> choice;

    if (choice == 'y' || choice == 'Y') {
        goto start;
    }

    cout << "Go read python books" << endl;
    return 0;
}
