#include <iostream>
#include <vector>
#include <math.h>
#pragma once

auto dot(const Field1D&a , const Field1D& b) -> double{
    double val = 0;
    for(int i = 0; i < a.size(); ++i){
        val+= (a.v[i]*b.v[i]);
    }
    return val;
}
auto axpy(const double a, const Field1D& x, Field1D& y) -> void{
    for(int i = 0; i < y.size(); ++i){
        y.v[i]+= (a*x.v[i]);
    }
}

auto norm2(const Field1D& x) -> double{ return sqrt(dot(x, x));}

// auto copy(const Field1D&, Field1D&) // TODO src, dst