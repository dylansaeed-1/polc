#include <cmath>
#include <functional>
#include <iostream>
#include <geometry.h>
#include <operators.h>
#include <numeric>
#include <linalg.h>
#include <solver.h>
#include <utils.h>
#include <test.h>
constexpr double PI = 3.14159265;
namespace Test {
    using Function =std::function<double(double x, double y)>;

void run_test_CG(int N){
    const Grid2D grid(N, N, 1.0, 1.0);
    Field2D K(grid, 1.0), q(grid, 0.0), p(grid, 0.0), g(grid, 0.0);

    // const Function p_sol = [&](double x, double y){
    //     return sin(PI*x)*sin(PI*y);
    // };
    Field2D p_exact(grid);
    for(int i = 0; i < grid.Nx; i++){
        double x = grid.x_cell(i);
        for(int j = 0; j < grid.Ny; j++){
            double y = grid.y_cell(j);
            p_exact(i,j) = sin(PI*x)*sin(PI*y);
            q(i, j) = 2* PI * PI *sin(PI*x)*sin(PI * y);
        }
    }

    const BC2D bc = set_bc(grid, p_exact);
    // const auto rhs = [&](double x, double y){return 2* PI * PI *sin(PI*x)*sin(PI * y);};
    // set_field(rhs, q);    
    
    build_bc_contrib(K, p, bc, g);
    const auto A = [&](const Field2D& x, Field2D& Ax){
        apply_A_2D(K, x, Ax);  // homogeneous operator (no g => SPD)
    };
    //Solving Ax = q-g
    int iters = cg_solve_2D(A, q-g, p, 5000, 1e-10);
    write_field_csv(K, "k_case1.csv");
    write_field_csv(p, "p_case1.csv");
    std::cout << "CG 2D converged in " << iters << " iterations\n";
}

void run_test_JFNK(int N){    
    const Grid2D grid(N, N, 1.0, 1.0);
    Field2D p(grid), q_nl(grid), p_exact(grid), R(grid), q_zero(grid);
    
    const Function sol_func = [&](double x, double y){return sin(PI*x)*sin(PI * y);};
    set_field(sol_func, p_exact);

    const auto bc = set_bc(grid, p_exact);
    richards_residual(p_exact, q_zero, bc, R);

    q_nl = R;
    
    auto F = [&](const Field2D& p, Field2D& r){ //Residual function we want to minimise
        richards_residual(p, q_nl, bc, r);
    };
    //Solving Ax = q-g
    int iters = jfnk_solve_2D(F, p);
    write_field_csv(p, "p_case1_jfnk.csv");
    std::cout << "N=" << N << "\n";
    std::cout << "JFNK 2D converged in " << iters << " Newton iterations\n";
    auto err = p - p_exact;
    std::cout << "L2 Error: " << norm2(err)*sqrt(grid.dx*grid.dy) << " " << grid.dx << "\n";
    
}

void run_test_richards(int N){
    
    int Nx = N, Ny = N;
    const Grid2D grid(Nx, Ny, 1.0, 1.0);
    Field2D p(grid, -0.5), q(grid, 0.0), R(grid);
    BC2D bc{
        BCSide2D::Dirichlet(Ny, 0.0),  // left
        BCSide2D::Neumann(Ny, 0.0),  // right
        BCSide2D::Neumann(Nx, 0.0),  // bottom
        BCSide2D::Dirichlet(Nx, -5.0)   // top
    };
    
    richards_residual(p, q, bc, R);
    auto F = [&](const Field2D& p, Field2D& r){ //Residual function we want to minimise
        richards_residual(p, q, bc, r);
    };
    int iters = jfnk_solve_2D(F, p);
    auto K = K_richards(p);
    write_field_csv(p, "p_richards.csv");
    write_field_csv(K, "k_richards.csv");
    F(p, R);
    std::cout << "N=" << N << "\n" << "JFNK 2D converged in " << iters << " Newton iterations\n";
    std::cout << "L2 Error: " << norm2(R)*sqrt(grid.dx*grid.dy) << " " << grid.dx << "\n";

}
}