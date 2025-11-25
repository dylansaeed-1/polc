#include "linalg.h"
#include "utils.h"
#include <functional>
#include <iostream>
#include <vector>
#include <cassert>
#include <geometry.h>

#pragma once

/*
Flux sign > 0 means flow to the +x direction

Darcy flux at face F = K_f* (p_r - p_l)/(dx)
p_r and p_l are adjacent cells or ghost pressures

K_{i+1/2} = 2*K_iK_i+1/(K_i*K_i+1) (harmonic mean)

operator value per cell is (Ap)_i = (F_{i+1/2} - F_{i-1/2}) /(dx)
N = 3 => | 0 | 1 | 2 |
*/

inline void apply_A_1D(const Field1D& K, const BC1D& bc, const Field1D& p, Field1D& Ap){
    //Sanity checks
    assert(K.size() == p.size());
    Ap.fill(0.0);
    
    const double dx = p.g.dx;
    const int n = K.g.N;

    //Compute left boundary cond.
    double f0 = 0.0;
    if (bc.left.type == BCSide1D::Type::Dirichlet){
        f0 = -K(0) * (2 * (p(0) - bc.left.val) / dx);   
    }
    else if (bc.left.type == BCSide1D::Type::Neumann){
        f0 = bc.left.val;
    }
    Ap(0)-= (f0/dx);

    //Compute Interior faces
    for (auto i = 0; i < K.size()-1; ++i){
        double F_r, K_r;
        K_r = 2.0*(K(i)*K(i+1))/(K(i) + K(i+1));
        F_r = -K_r * (p(i+1) - p(i)) / dx;

        Ap(i)+= (F_r / dx);
        Ap(i+1)-= (F_r /dx);

    }

    //Compute right boundary cond.
    double fn = 0.0;
    if (bc.right.type == BCSide1D::Type::Dirichlet){
        fn = -K(n-1) * (2 * (bc.right.val - p(n-1)) / dx);
    }
    else if (bc.right.type == BCSide1D::Type::Neumann){
        fn = bc.right.val;
    }

    Ap(n-1)+= (fn / dx);
    return;
}

// TODO(dylan): factor BC loops into a generic boundary helper once
// Darcy–Richards & Newton are in place.
inline void apply_A_2D(const Field2D& K, const Field2D& p, Field2D& Ap){
    assert(K.size() == p.size());
    Ap.fill(0.0);
 
    const double dx = p.g.dx, dy = p.g.dy;
    const int Nx = p.g.Nx, Ny = p.g.Ny;

    //Interior faces
    double F, K_r;
    for (int j = 0; j < Ny-1; ++j){
        for (int i = 0; i < Nx; ++i){            
            // //y flux
            K_r = 2.0 * K(i, j)*K(i, j+1) / (K(i, j) + K(i, j+1));
            F = -K_r * (p(i, j+1) - p(i, j))/dy;
            
            Ap(i, j)+= F/dy;
            Ap(i, j+1)-= F/dy;
        }
    }
    for (int i = 0; i < Nx - 1; ++i ){
        for (int j = 0; j < Ny; ++j){
            // x flux   
            K_r = 2.0 * K(i, j)*K(i+1, j) / (K(i, j) + K(i+1, j));
            F = -K_r * (p(i+1, j) - p(i, j))/dx;
            
            Ap(i, j)+= F/dx;
            Ap(i+1, j)-= F/dx;
        }
    }
        // Left boundary CURRENTLY ONLY HANDLES DIRICHLET
    for (int j = 0; j < Ny; ++j){
        K_r = K(0, j);
        F = -K_r * (2.0 * (p(0, j))/dx);
        Ap(0, j)-= F/dx;
    }

    //Bottom boundary (Dirichlet only for now)
    for (int i = 0; i < Nx; ++i){
        K_r = K(i, 0);
        F = -K_r * (2.0 * (p(i, 0))/dy);
        Ap(i, 0)-= F/dy;
    }

    // Right boundary (Dirichlet only for now)
    for (int j = 0; j < Ny; ++j){
        K_r = K(Nx-1, j);
        F = -K_r * (2.0 * (p(Nx-1, j))/dx);
        Ap(Nx-1, j)-= F/dx;
    }

    //Top boundary (Dirichlet only for now)
    for (int i = 0; i < Nx; ++i){
        K_r = K(i, Ny-1);
        F = -K_r * (2.0 * (p(i, Ny-1))/dy);
        Ap(i, Ny-1)-= F/dy;
    }
}
// g stores contribution of boundary flux to operator residual.
// Left/Bottom faces contribute with negative sign since flux leaves domain.
// Right/Top contribute with positive sign.
inline void build_bc_contrib(const Field2D& K,const Field2D& p, const BC2D& bc, Field2D& g){
    //Left boundary (Dirichlet only for now)
    double K_r, b_val, F;
    const int Nx = K.g.Nx, Ny = K.g.Ny;
    const double dx = K.g.dx, dy = K.g.dy;
    
    for (int j = 0; j < Ny; ++j){
        b_val = bc.left.val[j];
        if (bc.left.type == BCSide2D::Type::Dirichlet){
            K_r = K(0, j);
            F = -K_r * (2.0 * (-b_val)/dx);
            g(0, j)-= F/dx;
        }
        else{
            g(0,j)-=b_val/dx;
        }
    }
    //Bottom boundary (Dirichlet only for now)
    for (int i = 0; i < Nx; ++i){
        b_val = bc.bottom.val[i];
        if (bc.bottom.type == BCSide2D::Type::Dirichlet){
            K_r = K(i, 0);
            F = -K_r * (2.0 * (-b_val)/dy);
            g(i, 0)-= F/dy;
        }
        else{
            g(i, 0)-=b_val/dy;
        }
    }

    //Right boundary (Dirichlet only for now)
    for (int j = 0; j < Ny; ++j){
        b_val = bc.right.val[j];
        if (bc.right.type == BCSide2D::Type::Dirichlet){            
            K_r = K(Nx-1, j);
            F = -K_r * (2.0 * (b_val)/dx);
            g(Nx-1, j)+= F/dx;
        }
        else{
            g(Nx-1, j)+=b_val/dx;
        }
    }

    //Top boundary (Dirichlet only for now)
    for (int i = 0; i < Nx; ++i){
        b_val = bc.top.val[i];
        if (bc.top.type == BCSide2D::Type::Dirichlet){
            K_r = K(i, Ny-1);
            F = -K_r * (2.0 * (b_val)/dy);
            g(i, Ny-1)+= F/dy;
        }
        else{
            g(i, Ny-1)+= b_val/dy;
        }
    }
}

inline void jv_fd(
    const std::function<void(const Field2D& p, Field2D& F)> F, // Residual computation given current p iterate
    const Field2D& p,
    const Field2D& v,
    Field2D& Jv,
    Field2D& p_eps,
    Field2D& f0,
    Field2D& f_eps,    
    double eps=1e-6
){
    p_eps = p + eps*v;
    F(p, f0); 
    F(p_eps, f_eps);

    Jv = (f_eps - f0)/eps;
}


inline void richards_residual(const Field2D& p, const Field2D& q, const BC2D& bc, Field2D& R) {
    // 1. Compute K(p)  
    // 2. Compute Ap = ∇·(K(p) ∇p)
    // 3. R = Ap - q + g
    const auto& grid = p.g;
    const auto K = K_richards(p);
    
    Field2D Ap(grid, 0.0), g(grid, 0.0);
    apply_A_2D(K, p, Ap);
    build_bc_contrib(K, p, bc, g);
    R = Ap - q + g;
} 

inline void residual_2D(const Field2D& K, const Field2D& p, const Field2D& q, Field2D& r){
    apply_A_2D(K, p, r);
    r-=q;
}