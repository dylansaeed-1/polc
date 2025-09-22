#include <iostream>
#include <geometry.h>
#include <operators.h>
#include <math.h>
#include <numeric>
// #define PI 3.14159265

/*
Make residual_1D(...) (you’ve got it).

Tiny main():

grid [0,1], N=64, K=1, q=0

BC: Dirichlet p_L=0, p_R=1

init p_i = (i+0.5)/N

compute r = Ap − q

print ‖r‖∞ and ‖r‖₂

Expect ‖r‖∞ ~ 1e−12…1e−10.
*/
int main() {
    int lens[3] = {32, 64, 128};
    for (auto N : lens){
        auto g =  Grid1D(N, 1);
        const double PI = 3.141592653589793;
        auto bc = BC1D({BCSide::Dirichlet(0.0), BCSide::Dirichlet(0.0)});
    
        auto K = Field1D(g, 1.0);
        auto p = Field1D(g, 0.0);
        auto q = Field1D(g, 0.0);
        for (auto i = 0; i < p.v.size(); ++i){
            // p(i) = ((double)i + 0.5)/p.g.N; // Linear pressure field
            p(i) = sin(PI * g.x_cell(i));
            q(i) = (PI*PI)*sin(PI * g.x_cell(i));
        }
    
        double norm_inf = 0.0;
        auto r = Field1D(g, 0.0);
        
        residual_1D(K, p, bc, q, r);
        for (int i = 0; i < r.size(); ++i) {
            norm_inf = std::max(norm_inf, std::abs(r(i)));
        } 
        auto balance = std::accumulate(q.v.begin(), q.v.end(), 0.0);
    
        std::cout << "N " << N << ": " << norm_inf << std::endl;;
        // std::cout << "Bal. " << N << ": " << balance*p.g.dx << std::endl;;

    } 

    return 0;
}
