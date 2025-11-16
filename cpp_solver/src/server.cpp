#include <iostream>
#include <geometry.h>
#include <operators.h>
#include <math.h>
#include <numeric>
#include <linalg.h>
#include <cg.h>
#define PI 3.14159265

/*
Make residual_1D(...) (you’ve got it).

Tiny main():

grid [0,1], N=64, K=1, q=0

BC: Dirichlet p_L=0, p_R=1

init p_i = (i+0.5)/N

compute r = Ap − q  

print ‖r‖∞ and ‖r‖₂

Expect ‖r‖∞ ~ 1e−12…1e−10.
*/
int main() {

    int N = 50;
    double L = 1.0;
    Grid1D grid(N, L);

    // permeability K = 1   
    Field1D K(grid, 1.0);

    // RHS q = 0
    Field1D q(grid, 0.0);

    // initial guess p = 0
    Field1D p(grid, 0.0);


    BC1D bc{
        BCSide::Dirichlet(1.0), // left
        BCSide::Dirichlet(0.0)  // right
    };


    int max_iter = 1000;
    double tol = 1e-8;

    int iters = cg_solve_1D(K, bc, q, p, max_iter, tol);

    std::cout << "CG converged in " << iters << " iterations\n";

    // Compare with exact solution p(x) = 1 - x
    double max_err = 0.0;
    for (int i = 0; i < grid.N; ++i) {
        double x = grid.x_cell(i);
        double p_exact = 1.0 - x;
        double err = std::abs(p(i) - p_exact);
        max_err = std::max(max_err, err);
    }

    std::cout << "Max error = " << max_err << "\n";

    // for (auto N : lens){
    //     auto g =  Grid1D(N, 1);
    //     auto bc = BC1D({BCSide::Dirichlet(0.0), BCSide::Dirichlet(1.0)});
    
    //     auto K = Field1D(g, 1.0);
    //     auto p = Field1D(g, 0.0);
    //     auto q = Field1D(g, 0.0);
    //     for (auto i = 0; i < p.v.size(); ++i){
    //         // p(i) = ((double)i + 0.5)/p.g.N; // Linear pressure field
    //         // p(i) = sin(PI * g.x_cell(i));  //  Sinusoidal pressure field
    //         // q(i) = (PI*PI)*sin(PI * g.x_cell(i)); // true sol for p(i) above
    //         if (g.x_cell(i) < 0.5){K(i) = 1;}
    //         else {K(i) = 10;}
    //         p(i) = 1 - g.x_cell(i);
    //     }

    
    //     // cg_solve_1D(K, bc, q, p);
    //     auto r = Field1D(g, 0.0);
    //     residual_1D(K, bc, p, q, r);
    //     auto norm_inf = inf_norm(r);
        
    //     // auto balance = std::accumulate(q.v.begin(), q.v.end(), 0.0);
    
    //     std::cout << "N " << N << ": " << norm_inf << std::endl;;
    //     // std::cout << "Bal. " << N << ": " << balance*p.g.dx << std::endl;;

    // } 

    return 0;
}
