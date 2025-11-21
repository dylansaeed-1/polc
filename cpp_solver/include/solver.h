#include <iostream>
#include <vector>
#include <cassert>
#include <geometry.h>
#include <operators.h>
#include <limits>
#include <linalg.h>
#include <functional>
#pragma once

inline auto cg_solve_1D(
    const Field1D& K, 
    const BC1D& bc, 
    const Field1D& q, 
    Field1D& p, 
    int max_iter=100,
    double tol=1e-6
) -> int{

    int i = 0;
    auto grid = K.g;
    Field1D r(grid), z(grid), c(grid);

    apply_A_1D(K, bc, p, c);
    residual_1D(K, bc, p, q, r);
    scal(-1.0, r); //r = b - Ax
    
    auto d = r;
    double del_new = dot(r, r);
    double del_0 = del_new;
    while (i < max_iter and del_new > pow(tol, 2)*del_0){
        apply_A_1D(K, bc, d, z); // z = Ad
        axpy(-1.0,  c, z);
        auto alpha = del_new/dot(d, z);
        
        axpy(alpha, d, p); //x = x + alpha*d
        axpy(-1.0*alpha, z, r); //r = r - alpha*z

        auto del_old = del_new;
        del_new = dot(r, r);
        auto beta = del_new/del_old;

        scal(beta, d);
        axpy(1.0, r, d); // d = r + beta*d
        ++i;
    }
    return i;
}

inline auto cg_solve_2D(
    const std::function<void(const Field2D& x, Field2D& Ax)> A, // How to compute Ax
    const Field2D& b, // source term (Ax = b)
    Field2D& x,  // iniital guess
    int max_iter=100,
    double tol=1e-6
) -> int{

    int i = 0;
    auto grid = b.g;
    Field2D r(grid), z(grid), Ax(grid);

    A(x, Ax);
    r = b - Ax;    

    auto d = r;
    double del_new = dot(r, r);
    double del_0 = del_new;
    while (i < max_iter and del_new > pow(tol, 2)*del_0){
        A(d, z); // z = Ad
        auto alpha = del_new/dot(d, z);
        
        axpy(alpha, d, x); //x = x + alpha*d
        axpy(-1.0*alpha, z, r); //r = r - alpha*z

        auto del_old = del_new;
        del_new = dot(r, r);
        auto beta = del_new/del_old;

        scal(beta, d);
        axpy(1.0, r, d); // d = r + beta*d
        ++i;
    }
    return i;
}

inline int jfnk_solve_2D(
    std::function<void(const Field2D& p, Field2D& F)> F,
    Field2D& p,
    int max_newton = 20,
    double tol=1e-6
){
    Field2D R(p.g), dp(p.g), v(p.g ,0.0);
    /*
    Outer Newton loop
    Solves J(p_k)*v = -F(p_k)
    F(p_k) = residual at p_k
    J(p_k)*v is approximated using finite differences jp_v(p_k)
        Inner Krylov loop (CG)
        Solves J(p_k)*v = -F(p_k) using conjugate gradient method
    */
    for(int iter = 0; iter < max_newton; ++iter){
        F(p, R);
        double err = norm2(R);
        if (err < tol){
            return iter;
        }
        scal(-1.0, R);
        auto Jv = [&](const Field2D& v, Field2D& Jv_out){
            jv_fd(F, p, v, Jv_out);
        };
        dp.fill(0.0);
        cg_solve_2D(Jv, R, dp, 200, 1e-16);
        axpy(1.0, dp, p); // p_k+1 = p_k + dp
    }
    return max_newton;
}
