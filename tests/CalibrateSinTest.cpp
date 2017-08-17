#include "../include/ComputationalLib/CalibrateSin.hpp"
#include <iostream>

using namespace std;
using namespace ComputationalLib;
using namespace SimulationLib;

int main() {
    // BoundDouble x{RatioToNumeric<double>(PI_ratio{})};
    BoundDouble x{0};
    BoundDouble x_f{0};

    double standardDeviation = 1;
    StatisticalDistributions::Normal nDist(0, standardDeviation);

    // A Mersenne Twister pseudo-random generator of 64-bit 
    // numbers with a state size of 19937 bits.
    std::mt19937_64 rands(std::time(NULL));

    F f = [&] (BoundDouble x) {
        return std::sin(x()); // + nDist(rands); 
    }; 

    x_f = CalibrateSin(x, f);

    cout << x_f() << ":" << f(x_f()) << endl;
}
