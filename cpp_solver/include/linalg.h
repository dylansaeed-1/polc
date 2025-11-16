#include <iostream>
#include <vector>
#include <math.h>
#pragma once

//dot product
auto dot(const Field1D&a , const Field1D& b) -> double{
    assert(a.size() == b.size());
    
    double val = 0;
    for(int i = 0; i < a.size(); ++i){
        val+= (a(i)*b(i));
    }
    return val;
}

// y= y + ax
auto axpy(const double a, const Field1D& x, Field1D& y) -> void{
    assert(x.size() == y.size());
    
    for(int i = 0; i < y.size(); ++i){
        y(i)+= (a*x(i));
    }
}


auto norm2(const Field1D& x) -> double{ 
    return sqrt(dot(x, x));
}

// Src -> dst
auto copy(const Field1D& x, Field1D& y) -> void{
    assert(x.size() == y.size());
    
    for (int i = 0; i < x.size(); ++i){
        y(i) = x(i);
    }
} 

// Scaling  a * x
auto scal(const double a, Field1D& x) -> void{
    for (int i = 0; i < x.size(); ++i){
        x(i) *= a;
    }
}

auto inf_norm(Field1D& x) -> double{
    double norm_inf = 0.0;
    for(int i = 0; i < x.size(); ++i){
        norm_inf = std::max(norm_inf, x(i));
    }
    return norm_inf;
}