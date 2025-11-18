#include <geometry.h>
#include <linalg.h>
#include <operators.h>

#include <iostream>
#include <vector>

struct CGResult{
    int iters;              // iterations
    std::vector<double> sol; //Store solution
    std::vector<double> res_hist; // residual history
};

auto conjugate_gradient(const Field1D& K, const BC1D& bc, const Field1D& q, Field1D& p, int max_iter=100, double tol=1e-6) -> void{
    
}
