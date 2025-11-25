#include <iostream>
#include <vector>
#include <math.h>
#include <cassert>
#pragma once

//dot product
template <typename Field>
auto dot(const Field&a , const Field& b) -> double{
    assert(a.size() == b.size());
    
    double val = 0;
    for(int i = 0; i < a.size(); ++i){
        val+= (a.v[i]*b.v[i]);
    }
    return val;
}

// y= y + ax
template <typename Field>
auto axpy(const double a, const Field   & x, Field& y) -> void{
    assert(x.size() == y.size());
    
    for(int i = 0; i < y.size(); ++i){
        y.v[i]+= (a*x.v[i]);
    }
}

template <typename Field>
auto norm2(const Field& x) -> double{ 
    return sqrt(dot(x, x));
}

// Src -> dst
template <typename Field>
auto copy(const Field& x, Field& y) -> void{
    assert(x.size() == y.size());
    
    for (int i = 0; i < x.size(); ++i){
        y.v[i] = x.v[i];
    }
} 

// Scaling  a * x
template <typename Field>
auto scal(const double a, Field& x) -> void{
    for (auto& xi : x.v){
        xi *= a;
    }
}

template <typename Field>
auto inf_norm(const Field& x) -> double{
    double norm_inf = 0.0;
    for(auto xi : x.v){
        norm_inf = std::max(norm_inf, std::abs(xi));
    }
    return norm_inf;
}