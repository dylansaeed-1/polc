#include <iostream>
#include <vector>
#include <cassert>
#include <geometry.h>
#include <operators.h>
#include <limits>
#include <linalg.h>
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
    const Field2D& K, 
    // const BC2D& bc, 
    const Field2D& q, 
    Field2D& p, 
    int max_iter=100,
    double tol=1e-6
) -> int{

    int i = 0;
    auto grid = K.g;
    Field2D r(grid), z(grid), c(grid);

    // apply_A_2D(K, p, c);
    residual_2D(K, p, q, r);
    scal(-1.0, r); //r = b - Ax
    
    auto d = r;
    double del_new = dot(r, r);
    double del_0 = del_new;
    while (i < max_iter and del_new > pow(tol, 2)*del_0){
        apply_A_2D(K, d, z); // z = Ad
        // axpy(-1.0, c, z);
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