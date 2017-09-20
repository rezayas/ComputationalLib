#include "../include/ComputationalLib/CalibrateSinNDimensions.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace ComputationalLib;
using namespace SimulationLib;

#define OMEGA    0.001
#define ALPHA    0.05
#define EPSILON  0.00001

#define MAXITERS 1000
#define B_       50

template <typename... Xts,
          size_t...     I,
          typename     Xs = std::tuple<Xts...>,
          typename     Fn = std::function<double(Xts...)>,
Xs CalibrateSin(Xs x_i,
                Fn f,
                size_t d,
                std::index_sequence<I...>)
{
    auto Zero = IdentityFnGen<double, size_t>(0);

    Xs x_prev {Zero(I)...};

#ifdef DEBUG
    ofstream fout("CalibrateSin.csv");
    fout << "Iteration, Xs, f(Xs...)" << endl;
#endif

    double omega   {OMEGA};
    double alpha   {ALPHA};
    double epsilon {EPSILON};

    int maxIterations {MAXITERS};
    int b             {B_};
    int idx           {1};

    Xs lowers {std::get<I>(x_i).Lower...};
    Xs uppers {std::get<I>(x_i).Upper...};

    Polynomial g( Eigen::VectorXi::LinSpaced() )

}

template <typename... Xts,
          typename     Xs = std::tuple<Xts...>,
          typename     Fn = std::function<double(Xts...)>,
          size_t        d = std::tuple_size<Ns>,
          typename     IS = std::make_index_sequence<d> >
Xs CalibrateSin(Xs x_i, Fn f)
{
    return CalibrateSin(std::forward<Xs>(x_i),
                        std::forward<Fn>(f),
                        d
                        IS{});
}

constexpr std::tuple<Nts...>
VecToTpl (const std::vector<>)

template <typename...   Nts,
          size_t...     I,
          typename Ns = std::tuple<Nts...> >
std::vector<double>
TplToVec (const Ns boundNumbers,
          const std::index_sequence<I...>)
{
    return {std::get<I>(boundNumbers)()...}
}

template <typename... Nts,
          typename Ns = std::tuple<Nts...>
          size_t i    = std::tuple_size<Ns>,
          typename IS = std::make_index_sequence<i>>
std::vector<double>
TplToVec (const Ns boundNumbers)
{
    return TplToVec(std::forward<Ns>(boundNumbers),
                    IS{})
}

template <typename I,
          typename... Ts>
function<I(Ts...)>
IdentityFnGen(I value)
{
    return [=value] (Ts... params) -> I {
        return value;
    };
}