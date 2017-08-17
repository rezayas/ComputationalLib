#include "polyreg.hpp"

#include <cmath>

#include <Bound.h>
// need to include TimeSeries.h in order to use namespace SimulationLib
#include <TimeSeries.h>
#include <Normal.h>

using namespace std;
using namespace ComputationalLib;
using namespace SimulationLib;

using TWO_PI_ratio = std::ratio_divide<PI_ratio, std::ratio<1,2>>;

// for now BoundDouble is hard coded with bounds, 0, 2PI
using BoundDouble = Bound<double, std::ratio<0,1>, TWO_PI_ratio>;

using F = function<double(BoundDouble)>;

BoundDouble CalibrateSin(BoundDouble _x, F _f);