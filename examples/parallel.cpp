
#include "problems.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

int main()
{
    constexpr int nx = 255;
    constexpr int ny = 191;
    constexpr int repetitions = 3;
    constexpr double rtol = 1e-10;
    constexpr double atol = 1e-12;
    constexpr int max_iter = 200000;
    auto problem = make_problem(nx, ny, 2);
    const Grid a = problem.grid();
    const int N = nx * ny;
    std::vector<double> times;
    double baseline_time = 0.0;
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Grid: " << nx << " x " << ny
              << " (composite plate)\n";
    std::cout << "Tolerance: atol=" << atol
              << ", rtol=" << rtol << '\n';
    std::cout << "Repetitions: " << repetitions << '\n';
#ifndef _OPENMP
    std::cerr << "ERROR: OpenMP is not enabled. Recompile with -fopenmp.\n";
    return 1;
#else
    std::cout << "OpenMP enabled, version: " << _OPENMP << "\n";
#endif
    std::cout << "\nThreads,MedianSeconds,Iterations,Residual,"
                 "Converged,Speedup\n";
    for (int threads : {1, 2, 4}) {
        omp_set_dynamic(0);
        omp_set_num_threads(threads);
        times.clear();
        int iterations = 0;
        double residual = 0.0;
        bool converged = true;
        {
            std::vector<double> x(N, 0.0);
            Result r = solve_fast(
                a, problem.b.data(), x.data(),
                rtol, atol, max_iter);

            if (!r.converged) {
                std::cerr << "Warm-up solve failed to converge.\n";
                return 2;
            }
        }
        for (int rep = 0; rep < repetitions; ++rep) {
            std::vector<double> x(N, 0.0);
            const auto start = std::chrono::steady_clock::now();
            Result r = solve_fast(
                a, problem.b.data(), x.data(),
                rtol, atol, max_iter);
            const auto stop = std::chrono::steady_clock::now();
            times.push_back(
                std::chrono::duration<double>(stop - start).count());
            iterations = r.iterations;
            residual = r.residual;
            converged = converged && r.converged;
        }
        std::sort(times.begin(), times.end());
        const double median = times[times.size() / 2];
        if (!converged) {
            std::cerr << "Solver failed to converge with "
                      << threads << " threads.\n";
            return 3;
        }
        if (threads == 1)
            baseline_time = median;
        const double speedup = baseline_time / median;
        std::cout << threads << ','
                  << median << ','
                  << iterations << ','
                  << residual << ','
                  << (converged ? "yes" : "no") << ','
                  << speedup << '\n';
    }
    std::cout << "\nTiming includes solve_fast, including residual "
                 "checks; excludes grid construction, initialization, "
                 "and warm-up.\n";
    return 0;
}