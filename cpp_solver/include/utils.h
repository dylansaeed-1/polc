#include <cstdlib>
#include <fstream>
#include <functional>
#include <geometry.h>
#include <math.h>
#pragma once


inline void write_field_csv(const Field2D& f, const std::string& filename){
    auto out = std::ofstream(filename);
    for (int j = 0; j < f.g.Ny; ++j){
        for (int i = 0; i < f.g.Nx; ++i){
            out << f(i, j);
            if (i < f.g.Nx - 1){ out << ',';} //comma separate values
        }
        out << "\n";
    }
}

inline BC2D set_bc(const Grid2D& p, const std::function<double(int x, int y)> p_sol){
    const auto Nx = p.Nx;
    const auto Ny = p.Ny;
    BC2D bc{
        BCSide2D::Dirichlet(Ny, 0.0),  // left
        BCSide2D::Dirichlet(Ny, 0.0),  // right
        BCSide2D::Dirichlet(Nx, 0.0),  // bottom
        BCSide2D::Dirichlet(Nx, 0.0)   // top
    };
    for (int i = 0; i < Nx; ++i){
        bc.bottom.val[i] = p_sol(i, 0);
        bc.top.val[i] = p_sol(i, Ny-1);
        bc.left.val[i] = p_sol(0, i);
        bc.right.val[i] = p_sol(Nx-1, i);
    }
    return bc;
}

// van Genuchten model
inline double Sc(const double p, const double alpha=3.0, const double n=3.0){
    const double m  = 1.0 - 1.0/n;
    const double inner = pow(1.0 + alpha*abs(p), n);
    return pow(inner, -m);
}
inline double K_rel(const double p){
    const auto n = 2.0;
    const auto m = 1.0 - 1.0/n;
    const double inner = pow(1 - pow(Sc(p), 1.0/m), m);
    return pow(Sc(p), 0.5) * pow(inner, 2);
}

inline Field2D K_richards(
    const Field2D& P,
    const double Ks = 1e-4
){
    Field2D K(P.g);
    for(int i = 0; i < K.size(); ++i){
        double p =  P.v[i];
        const auto scp = Sc(p);
        K.v[i] = Ks * K_rel(scp);
    }
    return K;
}

inline void set_field(
    const std::function<double(double x, double y)>& F,
    Field2D& field
){
    for (int i = 0; i < field.g.Nx; ++i){
        const double x = field.g.x_cell(i);
        for (int j = 0; j < field.g.Ny; ++j){
            const double y = field.g.y_cell(j);
            field(i, j) = F(x, y); 
        }
    }
}