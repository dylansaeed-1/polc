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
void run_test_lienar(int N){
    
    Grid1D grid(N, 1.0);
    Field1D K(grid, 1.0), q(grid, 0.0), p(grid, 0.0);
    BC1D bc{ 
        BCSide1D::Dirichlet(1.0), 
        BCSide1D::Dirichlet(0.0)
    };

    int max_iter = 1000;
    double tol = 1e-8;

    int iters = cg_solve_1D(K, bc, q, p, max_iter, tol);
    
    double max_err = 0.0;
    for (int i = 0; i < grid.N; ++i) {
        double x = grid.x_cell(i);
        double p_exact = 1.0 - x;
        double err = std::abs(p(i) - p_exact);
        max_err = std::max(max_err, err);
    }
    
    std::cout << "Converged in : " << iters << " steps," << " Max error = " << max_err << "\n";
    // return iters;
}


int main() {
    
    int N_test[3] = {50, 100, 200};
    for (auto val : N_test){
        std::cout << val << " ";
        run_test_lienar(val);
    }   

    return 0;
}
