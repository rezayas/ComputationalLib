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
#define OMEGA    0.001
#define ALPHA    0.05
#define EPSILON  0.00001

#define B_       50
#define MAXITERS 1000

#include <cmath>

auto square = [] (auto n) -> decltype(auto) {
	return std::pow(n, 2);
};

auto cube = [] (auto n) -> decltype(auto) {
	return std::pow(n, 3);
};

// The following commented code is currently BROKEN:
//
// template <typename NumT>
// using CalibrateResults = std::pair<std::vector<NumT>, NumT>;

// template <typename... Xts,
// 		  size_t...   I,
// 		  typename    RT = double,
// 		  typename    Xs = std::tuple<Xts...>,
// 		  typename    F  = std::function<RT(Xts...)> >
// CalibrateResults<RT>
// TestND(Xs XInit, F f, std::index_sequence<I...>)
// {
// 	auto XFinal = PolyRegCal::forward<Xs>(XInit),
// 								std::forward<F>(f));
// 	auto YFinal = f(get<I>(XFinal)...);

// 	return {TplToVec(XFinal), YFinal};
// }

// template <typename... Xts,
// 		  typename    Xs = std::tuple<Xts...>,
// 		  typename    RT = double,
// 		  typename    F  = std::function<RT(Xts...)> >
// CalibrateResults<RT>
// TestND(Xs XInit, F f)
// {
// 	const auto ts = std::tuple_size<Xs>();

// 	return TestND(std::forward<Xs>(XInit),
// 				  std::forward<F>(f),
// 				  std::make_index_sequence<ts>{});
// }


TEST_CASE("ndimensional: f(x) = sin(x) + 1", "[PolyRegCal]") {
    auto X = std::make_tuple(UnitDouble{0.001});

    std::function<double(UnitDouble)> f = [&] (UnitDouble x) -> double {
        return std::sin( x() ) + 1;
    };

    auto XFinal = PolyRegCal(X, f, OMEGA, ALPHA, EPSILON, B_, MAXITERS);

    auto y_f = f(std::get<0>(XFinal));

    REQUIRE( y_f == Approx(2) );
}

using ScalarDouble = Bound<double, rZero, rOne>;

TEST_CASE("ndimensional: f(x, b) = sin(x) + b + 1", "[PolyRegCal]") {
	using Xs = std::tuple<UnitDouble, ScalarDouble>;
	using F = std::function<double(UnitDouble, ScalarDouble)>;

	Xs X {0.7, 0.7};

	F f = [&] (UnitDouble x, ScalarDouble b) -> double {
		return std::sin( x() ) + b() + 1;
	};

	auto XFinal = PolyRegCal(X, f, OMEGA, ALPHA, EPSILON, B_, MAXITERS);
	auto yFinal = f(get<0>(XFinal), get<1>(XFinal));

	REQUIRE( get<0>(XFinal)() == Approx(M_PI/2).margin(0.02) );
	REQUIRE( get<1>(XFinal)() == Approx(1) );

	REQUIRE( yFinal == Approx(3) );
	// REQUIRE( true );
}

TEST_CASE("ndimensional: sphere function n=2", "[PolyRegCal]") {
	double margin {0.0001};
	auto ApproxZero = Approx(0.).margin(margin);
	size_t n {2};

	using Domain = Bound<double, rInt<-2>, rInt<2>>;
	using Xs = std::tuple<Domain, Domain>;
	using F = std::function<double(Domain,Domain)>;

	Xs XInit {-2, -2};

	F f = [&] (auto a, auto b) {
		return -1 * (pow(a(), 2) + pow(b(), 2));
	};

	auto XFinal = PolyRegCal(XInit, f, OMEGA, ALPHA, EPSILON, B_, MAXITERS);
	auto YFinal = f(get<0>(XFinal), get<1>(XFinal));

	REQUIRE( get<0>(XFinal)() == ApproxZero );
	REQUIRE( get<1>(XFinal)() == ApproxZero );

	REQUIRE( YFinal == ApproxZero );
}

// TEST_CASE("ndimensional: sphere function with noise n=2", "[PolyRegCal]") {
// 	double margin {0.0001};
// 	auto ApproxZero = Approx(0.).margin(margin);

// 	std::random_device rd;
// 	std::mt19937 gen(rd());
// 	std::normal_distribution<> d(0, 0.1);

// 	using Domain = Bound<double, rInt<-2>, rInt<2>>;
// 	using Xs = std::tuple<Domain, Domain>;
// 	using F = std::function<double(Domain,Domain)>;

// 	Xs XInit {-2, -2};



// 	F f = [&] (auto a, auto b) {
// 		double noise {d(gen)};
// 		return noise + -1 * (pow(a(), 2) + pow(b(), 2));
// 	};

// 	auto XFinal = PolyRegCal(XInit, f);
// 	auto YFinal = f(get<0>(XFinal), get<1>(XFinal));

// 	REQUIRE( get<0>(XFinal)() == ApproxZero );
// 	REQUIRE( get<1>(XFinal)() == ApproxZero );

// 	REQUIRE( YFinal == ApproxZero );
// }

// TEST_CASE("ndimensional: Rosenbrock", "[PolyRegCal]") {
// 	double margin {0.0001};
// 	auto ApproxZero = Approx(0.).margin(margin);
// 	auto ApproxOne  = Approx(1.).margin(margin);
// 	size_t n {2};

// 	using Domain = Bound<double, rInt<-10>, rInt<10>>;
// 	using Xs = std::tuple<Domain, Domain>;
// 	using F = std::function<double(Domain,Domain)>;

// 	Xs XInit {0.95, 1.06};

// 	F f = [&] (auto _a, auto _b) {
// 		double a {_a()};
// 		double b {_b()};
// 		double prelim { 100*square(b-square(a)) + square(a-1) };

// 		return -1 * prelim;
// 	};

// 	auto XFinal = PolyRegCal, f);
// 	auto YFinal = f(get<0>(XFinal), get<1>(XFinal));

// 	REQUIRE( get<0>(XFinal)() == ApproxOne );
// 	REQUIRE( get<1>(XFinal)() == ApproxOne );

// 	REQUIRE( YFinal == ApproxZero );
// }

// TEST_CASE("ndimensional: beale's function", "[PolyRegCal]") {
// 	double margin {0.0001};

// 	auto ApproxN = [&] (double N) { return Approx(N).margin(margin); };

// 	using Domain = Bound<double, std::ratio<-9,2>, std::ratio<9,2>>;

// 	using Xs = std::tuple<Domain, Domain>;
// 	using F = std::function<double(Domain,Domain)>;

// 	Xs XInit {3, 0.49};

// 	F f = [&] (auto _x, auto _y) {
// 		double x {_x()};
// 		double y {_y()};

// 		double t1 {1.5-x+x*y};
// 		double t2 {2.25-x+x*square(y)};
// 		double t3 {2.625-x+x*cube(y)};

// 		return square(t1) + square(t2) + square(t3);
// 	};

// 	auto XFinal = PolyRegCal, f);
// 	auto YFinal = f(get<0>(XFinal), get<1>(XFinal));

// 	REQUIRE( get<0>(XFinal)() == ApproxN(3) );
// 	REQUIRE( get<1>(XFinal)() == ApproxN(3) );

// 	REQUIRE( YFinal == ApproxN(0.) );
// }