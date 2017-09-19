#include "polyreg.hpp"

#include <cmath>

#include <Bound.h>
// need to include TimeSeries.h in order to use namespace SimulationLib
#include <TimeSeries.h>
#include <Normal.h>

using namespace std;
using namespace ComputationalLib;
using namespace SimulationLib;

using rationalTwoPi = std::ratio_divide<PI_ratio, std::ratio<1,2>>;
using rationalZero  = std::ratio<0,1>;

// for now UnitDouble is hard coded with bounds, 0, 2PI
using UnitDouble = Bound<double,
                         rationalZero,
                         rationalTwoPi>;

using F = function<double(UnitDouble)>;

UnitDouble CalibrateSin(UnitDouble _x, F _f);