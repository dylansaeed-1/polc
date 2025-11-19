#include <iostream>
#include <geometry.h>
#include <operators.h>
#include <math.h>
#include <numeric>
#include <linalg.h>
#include <cg.h>
#include <utils.h>
#define PI 3.14159265

// void run_test_lienar(int N){
    
//     Grid1D grid(N, 1.0);
//     Field1D K(grid, 1.0), q(grid, 0.0), p(grid, 0.0);
//     BC1D bc{ 
//         BCSide1D::Dirichlet(1.0), 
//         BCSide1D::Dirichlet(0.0)
//     };

//     int max_iter = 1000;
//     double tol = 1e-8;

//     int iters = cg_solve_1D(K, bc, q, p, max_iter, tol);
    
//     double max_err = 0.0;
//     for (int i = 0; i < grid.N; ++i) {
//         double x = grid.x_cell(i);
//         double p_exact = 1.0 - x;
//         double err = std::abs(p(i) - p_exact);
//         max_err = std::max(max_err, err);
//     }
    
//     std::cout << "Converged in : " << iters << " steps," << " Max error = " << max_err << "\n";
//     // return iters;
// }

/*
2D test for manufactured solution of p(x, y) = x^2 + y^2
=> darcy operator K d/dx (p) = (2x, 2y) and so q(x,y) = -4.


*/
void run_test_2D(int N){
    
    int Nx = N, Ny = N;
    Grid2D grid(Nx, Ny, 1.0, 1.0);
    Field2D K(grid, 1.0), q(grid, -4.0), p(grid, 0.0), g(grid, 0.0);


    BC2D bc{
        BCSide2D::Dirichlet(Ny, 0.0),  // left
        BCSide2D::Dirichlet(Ny, 0.0),  // right
        BCSide2D::Dirichlet(Nx, 0.0),  // bottom
        BCSide2D::Dirichlet(Nx, 0.0)   // top
    };
    for (int j = 0; j < N; ++j) {
        double y = grid.y_face(j);
        bc.left.val[j]  = 0.0*0.0 + y*y;           // x=0
        bc.right.val[j] = 1.0*1.0 + y*y;           // x=1
    }
    for (int i = 0; i < N; ++i) {
        double x = grid.x_face(i);
        bc.bottom.val[i] = x*x + 0.0*0.0;          // y=0
        bc.top.val[i]    = x*x + 1.0*1.0;          // y=1
    }
    
    build_bc_contrib(K, bc, g);
    q = q - g;
    int iters = cg_solve_2D(K, q, p, 5000, 1e-10);
    write_field_csv(K, "k_case1.csv");
    write_field_csv(p, "p_case1.csv");
    std::cout << "CG 2D converged in " << iters << " iterations\n";
}

// void debug(){
//     int Nx = 32, Ny = 32;
//     Grid2D grid(Nx, Ny, 1.0, 1.0);
//     Field2D K(grid, 1.0), q(grid, 0.0), p(grid, 1.0), Ap(grid, 0.0);
    
//     BC2D bc{
//         BCSide2D::Dirichlet(Ny, 1.0),  // left
//         BCSide2D::Dirichlet(Ny, 1.0),  // right
//         BCSide2D::Dirichlet(Nx, 1.0),  // bottom
//         BCSide2D::Dirichlet(Nx, 1.0)   // top
//     };    
//     apply_A_2D(K, bc, p, Ap);
//     std::cout << inf_norm(Ap) << "\n";

// }
int main() {
    
    // int N_test[3] = {50, 100, 200};
    // for (auto val : N_test){
    //     std::cout << val << " ";
    //     run_test_lienar(val);
    // }   

    int N = 50;
    run_test_2D(N);

    // debug();
    return 0;
}
