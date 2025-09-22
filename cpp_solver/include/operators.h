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

void apply_A_1D(const Field1D& K, const Field1D& p, const BC1D& bc, Field1D& Ap){
    //Sanity checks
    assert(K.size() == p.size());
    Ap.fill(0.0);
    
    const double dx = p.g.dx;
    const int n = K.g.N;

    //Compute left boundary cond.
    double f0 = 0.0;
    if (bc.left.type == BCSide::Type::Dirichlet){
        f0 = -K(0) * (2 * (p(0) - bc.left.val) / dx);   
    }
    else if (bc.left.type == BCSide::Type::Neumann){
        f0 = -bc.left.val;
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
    if (bc.right.type == BCSide::Type::Dirichlet){
        fn = -K(n-1) * (2 * (bc.right.val - p(n-1)) / dx);
    }
    else if (bc.right.type == BCSide::Type::Neumann){
        fn = bc.right.val;
    }

    Ap(n-1)+= (fn / dx);
    return;
}

void residual_1D(const Field1D& K, const Field1D& p, const BC1D& bc, const Field1D& q, Field1D& r){
    apply_A_1D(K, p, bc, r);
    for(auto i = 0; i < r.v.size(); ++i){
        r(i)-= q(i);
    }
}