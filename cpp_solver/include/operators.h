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

void apply_A_1D(const Field1D& K, const BC1D& bc, const Field1D& p, Field1D& Ap){
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

void apply_A_2D(const Field2D& K, const BC2D& bc, const Field2D& p, Field2D& Ap){
    assert(K.size() == p.size());
    Ap.fill(0.0);
 
    const double dx = p.g.dx;
    const double dy = p.g.dy;

    const int Nx = p.g.Nx;
    const int Ny = p.g.Ny;

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
    
    //Left boundary (Dirichlet only for now)
    for (int j = 0; j < Ny; ++j){
        if (bc.left.type == BCSide2D::Type::Dirichlet){
            auto pL = bc.left.val[j];
            K_r = K(0, j);
            F = -K_r * (2.0 * (p(0, j) - pL)/dx);

            Ap(0, j)-= F/dx;
        }
    }

    //Bottom boundary (Dirichlet only for now)
    for (int i = 0; i < Nx; ++i){
        if (bc.bottom.type == BCSide2D::Type::Dirichlet){
            auto pB = bc.bottom.val[i];
            K_r = K(i, 0);
            F = -K_r * (2.0 * (p(i, 0) - pB)/dy);

            Ap(i, 0)-= F/dy;
        }
    }

    //Right boundary (Dirichlet only for now)
    for (int j = 0; j < Ny; ++j){
        if (bc.right.type == BCSide2D::Type::Dirichlet){
            auto pR = bc.right.val[j];
            K_r = K(Nx-1, j);
            F = -K_r * (2.0 * (p(Nx-1, j) - pR)/dx);

            Ap(Nx-1, j)+= F/dx;
        }
    }

    //Top boundary (Dirichlet only for now)
    for (int i = 0; i < Nx; ++i){
        if (bc.top.type == BCSide2D::Type::Dirichlet){
            auto pT = bc.top.val[i];
            K_r = K(i, Ny-1);
            F = -K_r * (2.0 * (p(i, Ny-1) - pT)/dy);

            Ap(i, Ny-1)+= F/dy;
        }
    }
}


void residual_1D(const Field1D& K, const BC1D& bc, const Field1D& p, const Field1D& q, Field1D& r){
    apply_A_1D(K, bc, p, r);
    for(auto i = 0; i < r.size(); ++i){
        r(i)-= q(i);
    }
}