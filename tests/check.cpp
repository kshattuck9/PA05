#include "sparse.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
void near(double x,double y) {
    if(!std::isfinite(x)||std::abs(x-y)>1e-9) throw std::runtime_error("incorrect value");
}
int main() {
    try {
        double d[]={4,4,4,4},w[]={0,-1,0,-1},e[]={-1,0,-1,0};
        double s[]={0,0,-1,-1},n[]={-1,-1,0,0};
        Grid a{2,2,d,w,e,s,n};
        double b[]={-1,3,7,11},x[]={1,2,3,4},y[4];
        apply(a,x,y); for(int p=0;p<4;++p) near(y[p],b[p]);
        double z[]={0,0,0,0};
        jacobi_sweep(a,b,z,y); for(int p=0;p<4;++p) near(y[p],b[p]/4);
        sor_sweep(a,b,z,1); near(z[0],-.25); near(z[1],.6875);
        near(z[2],1.6875); near(z[3],3.34375);
        for(double& v:z) v=0;
        auto r=solve(a,b,z,1,1.2,1e-12,1e-13,1000);
        if(!r.converged) throw std::runtime_error("did not converge");
        for(int p=0;p<4;++p) near(z[p],x[p]);
        // new test: fast solver on a known 2x2 system
        double x_fast[] = {0, 0, 0, 0};
        auto rf = solve_fast(a, b, x_fast, 1e-12, 1e-13, 1000);
        if (!rf.converged)
            throw std::runtime_error("fast solver did not converge");
        for (int p = 0; p < 4; ++p)
            near(x_fast[p], x[p]);
        double ax_fast[4];
        apply(a, x_fast, ax_fast);
        double residual_fast = 0.0;
        for (int p = 0; p < 4; ++p)
            residual_fast = std::max(
                residual_fast, std::abs(b[p] - ax_fast[p]));
        double b_norm = 0.0;
        for (int p = 0; p < 4; ++p)
            b_norm = std::max(b_norm, std::abs(b[p]));
        double tolerance = 1e-13 + 1e-12 * b_norm;
        if (residual_fast > tolerance)
            throw std::runtime_error("fast solver residual too large");
        near(residual_fast, rf.residual);
        std::cout<<"Public stencil, sweep, and solve checks passed.\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
