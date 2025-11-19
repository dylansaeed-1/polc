#include <fstream>
#include <geometry.h>
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