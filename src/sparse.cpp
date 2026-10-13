#include "sparse.h"
#include <limits>
#include <algorithm>
#include <cmath>

void apply(const Grid& a, const double* x, double* y) {
    const int N = a.nx * a.ny;
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < N; ++p) {
        const int i = p % a.nx;
        const int j = p / a.nx;
        y[p] = a.d[p] * x[p];
        if (i > 0)
            y[p] += a.w[p] * x[p - 1];
        if (i < a.nx - 1)
            y[p] += a.e[p] * x[p + 1];
        if (j > 0)
            y[p] += a.s[p] * x[p - a.nx];
        if (j < a.ny - 1)
            y[p] += a.n[p] * x[p + a.nx];
    } //TODO
}

void jacobi_sweep(const Grid& a, const double* b,
                  const double* old_x, double* new_x) {
    for (int p = 0; p < a.nx * a.ny; ++p) {
        int i = p % a.nx;
        int j = p / a.nx;
        double value = b[p];
        if (i > 0)
            value -= a.w[p] * old_x[p - 1];
        if (i < a.nx - 1)
            value -= a.e[p] * old_x[p + 1];
        if (j > 0)
            value -= a.s[p] * old_x[p - a.nx];
        if (j < a.ny - 1)
            value -= a.n[p] * old_x[p + a.nx];
        new_x[p] = value / a.d[p];
    }
 // TODO
}
void sor_sweep(const Grid& a, const double* b, double* x, double omega) {
    for (int p = 0; p < a.nx * a.ny; ++p) {
        int i = p % a.nx;
        int j = p / a.nx;
        double value = b[p];
        if (i > 0)
            value -= a.w[p] * x[p - 1];
        if (i < a.nx - 1)
            value -= a.e[p] * x[p + 1];
        if (j > 0)
            value -= a.s[p] * x[p - a.nx];
        if (j < a.ny - 1)
            value -= a.n[p] * x[p + a.nx];
        double x_gs = value / a.d[p];
        x[p] = (1.0 - omega) * x[p] + omega * x_gs;
    }
    // TODO
}
Result solve(const Grid& a, const double* b, double* x, int method, double omega, double rtol, double atol, int max_iter) {
    int N = a.nx * a.ny;
    double b_norm = 0.0;
    for (int p = 0; p < N; ++p)
        b_norm = std::max(b_norm, std::abs(b[p]));
    double* Ax = new double[N];
    double tolerance = atol + rtol * b_norm;
    apply(a, x, Ax);
    double residual = 0.0;
    for (int p = 0; p < N; ++p)
        residual = std::max(residual, std::abs(b[p] - Ax[p]));

    if (residual <= tolerance) {
        delete[] Ax;
        return {0, residual, true};
    }
    for (int iter = 1; iter <= max_iter; ++iter) {
        if (method == 0) {
            double* new_x = new double[N];
            jacobi_sweep(a, b, x, new_x);
            for (int p = 0; p < N; ++p)
                x[p] = new_x[p];

            delete[] new_x;
        }
        else {
            sor_sweep(a, b, x, omega);
        }
        apply(a, x, Ax);
        residual = 0.0;
        for (int p = 0; p < N; ++p)
            residual = std::max(residual, std::abs(b[p] - Ax[p]));
        if (residual <= tolerance) {
            delete[] Ax;
            return {iter, residual, true};
        }
    }
    delete[] Ax;
    return {max_iter, residual, false}; // TODO
}


Result solve_fast(const Grid& a, const double* b, double* x,
                  double rtol, double atol, int max_iter)
{
    const int N = a.nx * a.ny;
    double b_norm = 0.0;
    #pragma omp parallel for reduction(max:b_norm) schedule(static)
    for (int p = 0; p < N; ++p) {
        b_norm = std::max(b_norm, std::abs(b[p]));
    }
    const double tolerance = atol + rtol * b_norm;
    double* Ax = new double[N];
    apply(a, x, Ax);
    double residual = 0.0;
    #pragma omp parallel for reduction(max:residual) schedule(static)
    for (int p = 0; p < N; ++p) {
        residual = std::max(residual, std::abs(b[p] - Ax[p]));
    }
    if (residual <= tolerance) {
        delete[] Ax;
        return {0, residual, true};
    }
    const double pi = std::acos(-1.0);
    const double rho_j = 0.5 * (
        std::cos(pi / (a.nx + 1)) +
        std::cos(pi / (a.ny + 1))
    );
    const double omega =
        2.0 / (1.0 + std::sqrt(1.0 - rho_j * rho_j));
    constexpr int check_every = 5;
    int iterations = 0;
    bool converged = false;
    #pragma omp parallel shared(Ax, residual, iterations, converged)
    {
        for (int iter = 1; iter <= max_iter; ++iter) {
            #pragma omp for schedule(static)
            for (int p = 0; p < N; ++p) {
                const int i = p % a.nx;
                const int j = p / a.nx;
                if ((i + j) % 2 != 0)
                    continue;
                double value = b[p];
                if (i > 0)
                    value -= a.w[p] * x[p - 1];
                if (i + 1 < a.nx)
                    value -= a.e[p] * x[p + 1];
                if (j > 0)
                    value -= a.s[p] * x[p - a.nx];
                if (j + 1 < a.ny)
                    value -= a.n[p] * x[p + a.nx];
                const double x_gs = value / a.d[p];
                x[p] = (1.0 - omega) * x[p] + omega * x_gs;
            }
            #pragma omp for schedule(static)
            for (int p = 0; p < N; ++p) {
                const int i = p % a.nx;
                const int j = p / a.nx;
                if ((i + j) % 2 != 1)
                    continue;
                double value = b[p];
                if (i > 0)
                    value -= a.w[p] * x[p - 1];
                if (i + 1 < a.nx)
                    value -= a.e[p] * x[p + 1];
                if (j > 0)
                    value -= a.s[p] * x[p - a.nx];
                if (j + 1 < a.ny)
                    value -= a.n[p] * x[p + a.nx];
                const double x_gs = value / a.d[p];
                x[p] = (1.0 - omega) * x[p] + omega * x_gs;
            }
            if (iter % check_every == 0 || iter == max_iter) {
                #pragma omp single
                {
                    residual = 0.0;
                }
                #pragma omp for reduction(max:residual) schedule(static)
                for (int p = 0; p < N; ++p) {
                    const int i = p % a.nx;
                    const int j = p / a.nx;

                    double ap = a.d[p] * x[p];

                    if (i > 0)
                        ap += a.w[p] * x[p - 1];
                    if (i + 1 < a.nx)
                        ap += a.e[p] * x[p + 1];
                    if (j > 0)
                        ap += a.s[p] * x[p - a.nx];
                    if (j + 1 < a.ny)
                        ap += a.n[p] * x[p + a.nx];

                    Ax[p] = ap;
                    residual = std::max(residual, std::abs(b[p] - ap));
                }
                #pragma omp single
                {
                    iterations = iter;
                    converged = (residual <= tolerance);
                }
                if (converged)
                    break;
            }
        }
    }
    delete[] Ax;
    return {iterations, residual, converged};
}

/* Result solve_fast(const Grid& a, const double* b, double* x,
                  double rtol, double atol, int max_iter) {
    int N = a.nx * a.ny;
    double b_norm = 0.0;
    for (int p = 0; p < N; ++p)
        b_norm = std::max(b_norm, std::abs(b[p]));
    double tolerance = atol + rtol * b_norm;
    double* Ax = new double[N];
    apply(a, x, Ax);
    double residual = 0.0;
    for (int p = 0; p < N; ++p)
        residual = std::max(residual, std::abs(b[p] - Ax[p]));
    if (residual <= tolerance) {
        delete[] Ax;
        return {0, residual, true};
    }
    const double pi = std::acos(-1.0);

    const double rho_j =
        0.5 * (
            std::cos(pi / (a.nx + 1)) +
            std::cos(pi / (a.ny + 1))
        );

    const double omega =
        2.0 / (1.0 + std::sqrt(1.0 - rho_j * rho_j));

    for (int iter = 1; iter <= max_iter; ++iter) {
        for (int p = 0; p < N; ++p) {
            int i = p % a.nx;
            int j = p / a.nx;
            double value = b[p];
            if (i > 0)
                value -= a.w[p] * x[p - 1];
            if (i < a.nx - 1)
                value -= a.e[p] * x[p + 1];
            if (j > 0)
                value -= a.s[p] * x[p - a.nx];
            if (j < a.ny - 1)
                value -= a.n[p] * x[p + a.nx];
            double x_gs = value / a.d[p];
            x[p] = (1.0 - omega) * x[p] + omega * x_gs;
        }
        apply(a, x, Ax);
        residual = 0.0;
        for (int p = 0; p < N; ++p)
            residual = std::max(residual, std::abs(b[p] - Ax[p]));
        if (residual <= tolerance) {
            delete[] Ax;
            return {iter, residual, true};
        }
    }
    delete[] Ax;
    return {max_iter, residual, false};
    // Working fallback once solve is implemented. Optimize after correctness.
} */
