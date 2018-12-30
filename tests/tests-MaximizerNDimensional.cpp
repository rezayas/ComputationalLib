#include <tuple>
#include <random>

#include <Bound.h>

#include "catch.hpp"
#include "../include/ComputationalLib/PolyRegCal.hpp"

using namespace ComputationalLib;

using std::get;

using rTwoPi = std::ratio_divide<PI_ratio, std::ratio<1,2>>;
using rZero  = std::ratio<0,1>;
using rOne   = std::ratio<1,1>;

template <int n>
using rInt   = std::ratio<n,1>;

// for now UnitDouble is hard coded with bounds, 0, 2PI
using UnitDouble = Bound<double,
                         rZero,
                         rTwoPi>;

#define _USE_MATH_DEFINES // for C++

// Constants for algorithm
#define OMEGA    0.005
#define ALPHA    0.5

#define B_       100
#define MAXITERS 100
#define EPSILON  0.00001

#include <cmath>

auto square = [] (auto n) -> decltype(auto) {
	return std::pow(n, 2);
};

auto cube = [] (auto n) -> decltype(auto) {
	return std::pow(n, 3);
};

TEST_CASE("ndimensional: sphere function n=2", "[PolyRegCalMaximizer]") {
	double margin {0.0001};
	auto ApproxX0 = Approx(-1.).margin(margin);
	auto ApproxX1 = Approx(0.).margin(margin);
	auto ApproxY = Approx(0.).margin(margin);

	size_t n {2};

	using Domain = Bound<double, rInt<-5>, rInt<5>>;
	using Xs = std::tuple<Domain, Domain>;
	using F = std::function<double(Domain,Domain)>;

	std::default_random_engine generator;
    std::normal_distribution<double> dist(0.0, 1.0);

	Xs XInit {3, 2};

	F f = [&] (auto a, auto b) {
		// double x2 = b();
		// double penalty;
		// if (x2 < 2) {
		// 	penalty = - 1000* pow(min(x2-2, 0.0),2); 
		// 	x2 = 2;
		// }

		return -1 * (pow(a()+1, 2) + pow(b(), 2)) + dist(generator); 
	};

	auto XFinal = PolyRegCal(XInit, f, OMEGA, ALPHA, EPSILON, B_, MAXITERS);
	auto YFinal = f(get<0>(XFinal), get<1>(XFinal));

	REQUIRE( get<0>(XFinal)() == ApproxX0 );
	REQUIRE( get<1>(XFinal)() == ApproxX1 );

	REQUIRE( YFinal == ApproxY );
}
