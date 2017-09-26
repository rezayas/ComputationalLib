#include <iostream>
#include <fstream>
#include <utility>
#include <tuple>

#include <Eigen/core>

#include "polyreg.hpp"

// Comment the following line to disable logging
#define DEBUG

// Constants for algorithm
#define OMEGA    0.001
#define ALPHA    0.05
#define EPSILON  0.00001

#define B_       50
#define MAXITERS 1000

namespace ComputationalLib {

    using std::get;
    using std::min;
    using std::max;
    using std::sqrt;
    using std::pow;
    using std::endl;

    template <typename... Xts,
              typename     Xs = std::tuple<Xts...>,
              size_t...     I>
    double
    SubtractedNorm(Xs a, Xs b, std::index_sequence<I...>)
    {
        auto subtracted =
            std::make_tuple(get<I>(a)() - get<I>(b)()...);

        auto norm = sqrt( (... + pow(get<I>(subtracted), 2)) );

        return norm;
    }

    template <typename... Xts,
              typename     Xs = std::tuple<Xts...> >
    double
    SubtractedNorm(Xs a, Xs b)
    {
        const auto ts = std::tuple_size<Xs>();

        return SubtractedNorm(std::forward<Xs>(a),
                              std::forward<Xs>(b),
                              std::make_index_sequence<ts>{});
    }

    template <typename VectorT = double,
              typename...  Xts,
              typename      Xs = std::tuple<Xts...>,
              size_t...      I >
    std::vector<VectorT>
    TplToVec(Xs tpl, std::index_sequence<I...>)
    {
        return {get<I>(tpl)()...};
    }

    template <typename VectorT = double,
              typename...  Xts,
              typename      Xs = std::tuple<Xts...> >
    std::vector<VectorT>
    TplToVec(Xs tpl)
    {
        const auto ts = std::tuple_size<Xs>();

        return TplToVec(std::forward<Xs>(tpl),
                        std::make_index_sequence<ts>{});
    }

    // Returns a function 'f: Ts... => I' whose output
    // is always 'value'
    template <typename I,
              typename... Ts>
    std::function<I(Ts...)>
    IdentityFnGen(I value)
    {
        return [=] (Ts&&... params) -> I {
            return value;
        };
    }

    template <typename... Xts,
              size_t...     I,
              typename     FT = double,
              typename     Xs = std::tuple<Xts...>,
              typename     Fn = std::function<FT(Xts...)> >
    Xs
    PolyRegCal(Xs x_i,
                 const Fn &f,
                 const size_t &d,
                 std::index_sequence<I...>)
    {
        // Create a function 'Zero: size_t => double' which always
        // returns zero
        auto Zero = IdentityFnGen<FT, size_t>(0);

        // Create a tuple of 'previous' values, which in this case
        // we assume to be zero
        Xs x_prev {Zero(I)...};

#ifdef DEBUG
        // Open a CSV file to record the progress of the calibration
        // algorithm
        std::ofstream fout("PolyRegCal.csv");
        fout << "Iteration, [Xs], f(Xs...)" << endl;
#endif

        // Omega: The penalty for regressing towards steeper slopes,
        //   used in the construction of the polynomial used in regression
        // Alpha: The step size. It is multiplied by the slope to calculate
        //   the distance to the next vector x_i
        // Epsilon: The minimum distance for each step. When a step becomes
        //   smaller than epsilon, the algorithm terminates
        // 'b': A constant used to calculate Lambda, the learning weight.
        //   As 'b' increases, lambda approaches 1. Therefore, the values of
        //   'b' and 'maxIters' should be chosen such that "b is very
        //   close to 1 as b approaches 'maxIters'??
        double omega   {OMEGA};
        double alpha   {ALPHA};
        double epsilon {EPSILON};

        int b          {B_};
        int maxIters   {MAXITERS};
        int idx        {1};

        // Form two tuples of lower and upper bounds on each x in X
        Xs lowers {get<I>(x_i).Lower...};
        Xs uppers {get<I>(x_i).Upper...};

        // Set up a regression using a 2nd-degree polynomial for a
        // function f(X) where the set X has 'd' members
        PolynomialRegression Regression(2, d, omega);

        // The algorithm must stop iterating when we hit 'maxIters', or
        // when the norm of '(x_i - x_prev)' is smaller than 'epsilon'
        auto terminationPred = [&] (const auto &idx,
                                    const auto &x_i,
                                    const auto &x_prev,
                                    const auto &epsilon) -> bool {

            bool hitMaxIters {idx >= maxIters};
            bool insignificantDifference { SubtractedNorm(x_i, x_prev) < epsilon };

            // Terminate if either of these conditions is satisfied
            return hitMaxIters || insignificantDifference;
        };

        do {
            // Calculate learning weight for regression
            double lambda_i { (double)idx / (double)(b + idx) };
            double y_i      { f(get<I>(x_i)...) };

            // Represent x_i as an Eigen vector for use with Regression obj.
            // To do this, we create a temporary Eigen vector and assign the
            // contents of 'x_i' to it.
            auto x_i_V = TplToVec(x_i);
            Eigen::VectorXd x_i_E(d);
            for (size_t i = 0; i < x_i_V.size(); ++i)
                x_i_E[i] = x_i_V[i];
            // for (auto i = x_i_V.begin(); i != x_i_V.end(); ++i)
                // x_i_E << *i;


            // Add a new point to the regression
            Regression.updateCoefficients(x_i_E, y_i, lambda_i);

            // Numerically calculate the derivative of each value x in X at the
            // point X
            Eigen::VectorXd slope_i_E = Regression.poly().derivative(x_i_E);

            // Prepare for modification of 'x_i'
            x_prev = x_i;

            // For each x in X, travel along the slope a little bit in
            // the positive direction
            x_i = { max(get<I>(lowers)(),
                        min(get<I>(uppers)(),
                            get<I>(x_i)() + (alpha*slope_i_E[I])) )... };

#ifdef DEBUG
            fout << idx;
            (fout << ... << ("," + std::to_string(get<I>(x_prev)())));
            fout << ",ENDX," << y_i << endl;
#endif

            idx += 1;
        } while (!terminationPred(idx, x_i, x_prev, epsilon));

#ifdef DEBUG
        fout.flush();
        fout.close();
#endif

        return x_i;
    }

    template <typename... Xts,
              typename     FT = double,
              typename     Xs = std::tuple<Xts...> >
    Xs PolyRegCal(const Xs &x_i, const std::function<FT(Xts...)> &f)
    {
        const auto ts = std::tuple_size<Xs>();

        return PolyRegCal(std::forward<decltype(x_i)>(x_i),
                             std::forward<decltype(f)>(f),
                             ts,
                             std::make_index_sequence<ts>{});
    }

}
// Clean up constants
#undef OMEGA
#undef ALPHA
#undef EPSILON
#undef B_
#undef MAXITERS