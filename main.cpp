#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

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
            cerr << "Error: the function does not change the sign to [" << a << ", " << b << "]" << endl;
            iterations = 0;
            return NAN;
        }

        double c;
        iterations = 0;

        while ((b - a) / 2.0 > eps && iterations < max_iter) {
            c = (a + b) / 2.0;

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

start: // ← метка для возврата

    cout << "Enter the left border  a: ";
    cin >> left;
    cout << "Enter the right border b: ";
    cin >> right;
    cout << "Enter eps accuracy: ";
    cin >> eps;

    solver.setInterval(left, right);
    solver.setEpsilon(eps);
    solver.setMaxIter(1000);

    double root = solver.solve();

    if (!isnan(root)) {
        cout << fixed << setprecision(6);
        cout << "Root found: x = " << root << endl;
        cout << "f(x) = " << (root*root*root - 4*root*root + 5*root - 3) << endl;
        cout << "Iterations: " << solver.getIterations() << endl;
    }

    cout << "\nWould you like to solve another equation? (y/n): ";
    cin >> choice;

    if (choice == 'y' || choice == 'Y') {
        goto start;
    }

    cout << "Go read python books" << endl;
    return 0;
}
