#include "../include/ComputationalLib/CalibrateSin.hpp"
#include <iostream>

using namespace std;
using namespace ComputationalLib;
using namespace SimulationLib;

int main(int argc, char const *argv[]) {
    if (argc < 2) {
        cout << "Please give a starting value for x (between 0 and 2π)"
             << endl;
        exit(1);
    }


    // x is the value that will be given to CalibrateSin function
    UnitDouble x{stod(argv[1])};

    // x_f will be set to the final outputted value of the model
    UnitDouble x_f{0};

    // To add some noise, you can add values from a normal distribution
    double standardDeviation = 1;
    StatisticalDistributions::Normal nDist(0, standardDeviation);

    // A Mersenne Twister pseudo-random generator of 64-bit
    // numbers with a state size of 19937 bits.
    std::mt19937_64 rands(std::time(NULL));

    F f = [&] (UnitDouble x) {
        return std::sin(x()) + 1; // + nDist(rands);
    };

    x_f = CalibrateSin(x, f);

    cout << "Final x: "<< x_f() << endl;
    cout << "Final y: " << f(x_f()) << endl;
}
