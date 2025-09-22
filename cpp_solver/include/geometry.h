#include <iostream>
#include <vector>

#pragma once


// Struct for storing grid data
struct Grid1D{
    int     N;  // Number of cells
    double  L;  // domain length [0, L]
    double  dx; // spacing

    Grid1D(int N_, double L_) : N(N_), L(L_), dx(L_/N_){}
    inline double x_cell(int i) const{
        return (i + 0.5) * dx;
    }
    inline double x_face(int f) const{
        return (f * dx);
    }

};


//Storing Scalar fields
struct Field1D{
    const Grid1D& g;
    std::vector<double> v;
    explicit Field1D(Grid1D& grid, double init=0.0) : g(grid), v(static_cast<size_t>(grid.N), init){}

    double& operator()(int i){ return v[static_cast<size_t>(i)];}
    double operator()(int i) const {return v[static_cast<size_t>(i)];}
    inline int size() const {return g.N;}
    inline void fill(const double val){std::fill(v.begin(), v.end(), val);}
};

//Storing boundary conditions
struct BCSide{
    enum class Type {Dirichlet, Neumann};
    Type type;
    double val; //pressure if Dirichlet, flux if Neumann

    static BCSide Dirichlet(double p){
        return {Type::Dirichlet, p};
    }
    static BCSide Neumann(double g){
        return {Type::Neumann, g};
    }
};

struct BC1D{ //  boundary condiiton object
    BCSide left;
    BCSide right;
};
