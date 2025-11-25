#include <chrono>
#include <test.h>
#include <iostream>

int main() {    
    for(auto v : {16, 32, 64}){
        const auto start_time = std::chrono::high_resolution_clock::now();
        Test::run_test_richards(v);
        const auto end_time = std::chrono::high_resolution_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Solved in: " << duration.count() << "ms\n";
    }
    return 0;
}
