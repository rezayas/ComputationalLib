#include <tuple>

#include <Bound.h>

#include "catch.hpp"
#include "../include/ComputationalLib/CalibrateSinNDimensional.hpp"

using namespace ComputationalLib;

using std::get;

using rationalTwoPi = std::ratio_divide<PI_ratio, std::ratio<1,2>>;
using rationalZero  = std::ratio<0,1>;

// for now UnitDouble is hard coded with bounds, 0, 2PI
using UnitDouble = Bound<double,
                         rationalZero,
                         rationalTwoPi>;

TEST_CASE("Hello_world_CalibrateSinNDimensional", "[calibrateSinNDimensional]") {
    auto X = std::make_tuple(UnitDouble{0});

    std::function<double(UnitDouble)> f = [&] (UnitDouble x) -> double {
        return std::sin( x() ) + 1;
    };

    auto xFinal = CalibrateSinN(X, f);

    auto y_f = f(std::get<0>(xFinal));

    REQUIRE( y_f == Approx(2) );
}