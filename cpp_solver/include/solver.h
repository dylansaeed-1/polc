#include <iostream>
#include <vector>
#include <cassert>
#include <geometry.h>
#include <operators.h>
#include <limits>
#include <linalg.h>
#include <functional>
#pragma once

inline auto cg_solve_2D(
    const std::function<void(const Field2D& x, Field2D& Ax)> A, // How to compute Ax
    const Field2D& b, // source term (Ax = b)
    Field2D& x,  // iniital guess
    int max_iter=100,
    double tol=1e-6
) -> int{

    int i = 0;
    const auto grid = b.g;
    Field2D r(grid), z(grid), Ax(grid);

    A(x, Ax);
    r = b - Ax;    

    auto d = r;
    double del_new = dot(r, r);
    double del_0 = del_new;
    while (i < max_iter and del_new > tol*tol*del_0){
        A(d, z); // z = Ad
        const auto alpha = del_new/dot(d, z);
        
        x+=alpha*d;
        r-=alpha*z;
        const auto del_old = del_new;
        del_new = dot(r, r);
        const auto beta = del_new/del_old;

        d = r + beta*d;
        ++i;
    }
    return i;
}

/*
Outer Newton loop
Solves J(p_k)*v = -F(p_k)
F(p_k) = residual at p_k
J(p_k)*v is approximated using finite differences jp_v(p_k)
    Inner Krylov loop (CG)
    Solves J(p_k)*v = -F(p_k) using conjugate gradient method


TODO: Implement inexact Newton
*/
inline int jfnk_solve_2D(
    const std::function<void(const Field2D& p, Field2D& F)>& F,
    Field2D& p,
    int max_newton = 20,
    double tol=1e-6
){
    Field2D R(p.g), dp(p.g), v(p.g), p_eps(p.g), f0(p.g), f_eps(p.g);
    for(int iter = 0; iter < max_newton; ++iter){
        F(p, R);
        double err = norm2(R);
        if (err < tol){
            return iter;
        }
        R*=-1;
        auto Jv = [&](const Field2D& v, Field2D& Jv_out){
            jv_fd(F, p, v, Jv_out, p_eps, f0, f_eps);

        };
        dp.fill(0.0);
        cg_solve_2D(Jv, R, dp, 200, 1e-6);
        p+=dp;
    }
    return max_newton;
}
