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

struct Grid2D{
    int Nx, Ny;
    double Lx, Ly;
    double dx, dy;

    //Maybe helper functions
};


//Storing Scalar fields
struct Field1D{
    const Grid1D& g;
    std::vector<double> v;
    explicit Field1D(Grid1D& grid, double init=0.0) : g(grid), v(static_cast<size_t>(grid.N), init){}
    Field1D(Field1D& x): g(x.g), v(x.v){}

    Field1D& operator=(Field1D other){
        std::swap(v, other.v);
        return *this;
    }
    double& operator()(int i){ return v[static_cast<size_t>(i)];}
    double operator()(int i) const {return v[static_cast<size_t>(i)];}
    inline int size() const {return g.N;}
    inline void fill(const double val){std::fill(v.begin(), v.end(), val);}
};

struct Field2D{
    const Grid2D g;
    std::vector<double> v;
    explicit Field2D(Grid2D& grid, double init=0.0) : g(grid), v(static_cast<size_t>(grid.Nx*grid.Ny), init) {};
    
    Field2D(Field2D& x) : g(x.g), v(x.v) {};
    Field2D& operator=(Field2D other){
        std::swap(v, other.v);
        return *this;
    }
    inline int idx(int i, int j) const {return g.Nx*j + i;} //i = x, j = y
    double& operator()(int i, int j){
        return v[this->idx(i, j)];
    }
    const double& operator()(const int i, const int j) const{
        return v[this->idx(i, j)];
    }
    inline int size() const {return (int)v.size();}
    inline void fill(const double val){ std::fill(v.begin(), v.end(), val);} 
};


//Storing boundary conditions
struct BCSide1D{
    enum class Type {Dirichlet, Neumann};
    Type type;
    double val; //pressure if Dirichlet, flux if Neumann

    static BCSide1D Dirichlet(double p){
        return {Type::Dirichlet, p};
    }
    static BCSide1D Neumann(double g){
        return {Type::Neumann, g};
    }
};

struct BCSide2D{
    enum class Type {Dirichlet, Neumann};
    Type type;
    std::vector<double> val; //pressure if Dirichlet, flux if Neumann

    static BCSide2D Dirichlet(std::vector<double>& p){
        return {Type::Dirichlet, p};
    }
    static BCSide2D Neumann(std::vector<double>& g){
        return {Type::Neumann, g};
    }
};

struct BC2D{
    BCSide2D left, right, bottom, top;
};

struct BC1D{ //  boundary condiiton object
    BCSide1D left, right;
};
