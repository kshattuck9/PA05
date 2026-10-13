
#include "problems.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

struct Benchmark {
    std::string plate;
    int nx, ny;
    std::string method;
    double omega;
    int iterations;
    double seconds;
    double residual;
    bool converged;
};

Benchmark run_case(int nx, int ny, int material,
                   const std::string& method, int method_id,
                   double omega)
{
    auto p = make_problem(nx, ny, material);
    const Grid a = p.grid();
    const int N = nx * ny;

    std::vector<double> x(N, 0.0);
    std::vector<double> Ax(N, 0.0);

    const auto start = std::chrono::steady_clock::now();

    Result r = solve(a, p.b.data(), x.data(),
                     method_id, omega,
                     1e-10, 1e-12, 200000);

    const auto stop = std::chrono::steady_clock::now();
    const double seconds =
        std::chrono::duration<double>(stop - start).count();

    // Independently recompute the infinity-norm residual.
    apply(a, x.data(), Ax.data());

    double residual = 0.0;
    for (int k = 0; k < N; ++k) {
        residual = std::max(
            residual, std::abs(p.b[k] - Ax[k]));
    }

    return {
        material == 1 ? "Uniform" : "Composite",
        nx, ny, method, omega, r.iterations,
        seconds, residual, r.converged
    };
}

void write_temperature(const std::string& filename,
                       int nx, int ny)
{
    auto p = make_problem(nx, ny, 2);
    const Grid a = p.grid();
    const int N = nx * ny;

    std::vector<double> x(N, 0.0);
    std::vector<double> Ax(N, 0.0);

    // Use SOR, omega = 1.5, for the composite plate.
    Result r = solve(a, p.b.data(), x.data(),
                     1, 1.5, 1e-10, 1e-12, 200000);

    apply(a, x.data(), Ax.data());

    double residual = 0.0;
    for (int k = 0; k < N; ++k)
        residual = std::max(
            residual, std::abs(p.b[k] - Ax[k]));

    if (!r.converged) {
        std::cerr << "Temperature solve did not converge: "
                  << "iterations=" << r.iterations
                  << ", residual=" << residual << '\n';
        return;
    }

    std::ofstream out(filename);
    out << "x,y,T\n" << std::setprecision(17);

    const double hx = 1.0 / (nx + 1);
    const double hy = 1.0 / (ny + 1);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int k = j * nx + i;
            out << (i + 1) * hx << ','
                << (j + 1) * hy << ','
                << x[k] << '\n';
        }
    }

    std::cout << "\nTemperature grid " << nx << "x" << ny
              << ": residual=" << residual
              << ", iterations=" << r.iterations << '\n';
    std::cout << "Wrote " << filename << '\n';
}

int main()
{
    std::vector<Benchmark> results;

    const int sizes[] = {31, 63};

    for (int material : {1, 2}) {
        for (int m : sizes) {
            // Jacobi
            results.push_back(run_case(
                m, m, material, "Jacobi", 0, 1.0));

            // Gauss-Seidel is SOR with omega = 1.
            results.push_back(run_case(
                m, m, material, "GS", 1, 1.0));

            // SOR with the requested parameters.
            results.push_back(run_case(
                m, m, material, "SOR", 1, 1.5));

            results.push_back(run_case(
                m, m, material, "SOR", 1, 1.8));
        }
    }

    std::cout << '\n'
              << std::left
              << std::setw(11) << "Plate"
              << std::setw(9) << "Grid"
              << std::setw(9) << "Method"
              << std::setw(9) << "Omega"
              << std::right
              << std::setw(12) << "Iterations"
              << std::setw(12) << "Seconds"
              << std::setw(14) << "Residual"
              << "  Status\n";

    std::cout << std::scientific << std::setprecision(3);

    for (const auto& r : results) {
        std::cout << std::left
                  << std::setw(11) << r.plate
                  << std::setw(9) << (std::to_string(r.nx)
                                     + "x" + std::to_string(r.ny))
                  << std::setw(9) << r.method
                  << std::setw(9) << r.omega
                  << std::right
                  << std::setw(12) << r.iterations
                  << std::setw(12) << r.seconds
                  << std::setw(14) << r.residual
                  << "  " << (r.converged ? "OK" : "FAILED")
                  << '\n';
    }

    write_temperature("build/composite_temperature.csv",
                      127, 95);

    return 0;
}