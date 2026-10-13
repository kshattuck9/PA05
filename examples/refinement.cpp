
#include "problems.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

int main()
{
    const int grids[] = {15, 31, 63, 127};
    const double pi = std::acos(-1.0);
    double previous_error = 0.0;
    std::cout << std::scientific << std::setprecision(8);
    std::cout << "m,h,residual,max_error,discrete_error,ratio,iterations\n";
    for (int m : grids) {
        const double h = 1.0 / (m + 1);
        auto p = make_problem(m, m, 0);
        const Grid a = p.grid();
        const int N = m * m;
        std::vector<double> x(N, 0.0);
        std::vector<double> Ax(N, 0.0);
        Result r = solve_fast(
            a, p.b.data(), x.data(),
            1e-11, 1e-13, 1000000
        );
        apply(a, x.data(), Ax.data());
        double residual = 0.0;
        double max_error = 0.0;
        double discrete_error = 0.0;
        const double lambda_h =
            8.0 / (h * h) *
            std::pow(std::sin(pi * h / 2.0), 2);
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < m; ++i) {
                const int k = j * m + i;
                const double xx = (i + 1) * h;
                const double yy = (j + 1) * h;
                const double u =
                    std::sin(pi * xx) * std::sin(pi * yy);
                const double U =
                    (2.0 * pi * pi / lambda_h) * u;

                residual = std::max(
                    residual, std::abs(p.b[k] - Ax[k])
                );

                max_error = std::max(
                    max_error, std::abs(x[k] - u)
                );

                discrete_error = std::max(
                    discrete_error, std::abs(x[k] - U)
                );
            }
        }
        std::cout << m << ','
                  << h << ','
                  << residual << ','
                  << max_error << ','
                  << discrete_error << ',';
        if (previous_error > 0.0)
            std::cout << previous_error / max_error;
        else
            std::cout << "N/A";
        std::cout << ',' << r.iterations << '\n';
        if (!r.converged) {
            std::cerr << "Solver failed to converge for m="
                      << m << ", residual=" << residual << '\n';
            return 1;
        }
        previous_error = max_error;
    }
    return 0;
}
