#pragma once
/*
 
 The quadprog_solve() function implements the algorithm of Goldfarb and Idnani 
 for the solution of a (convex) Quadratic Programming problem
 by means of an active-set dual method.
	 
The problem is in the form:

min 0.5 * x G x + g0 x
s.t.
    CE^T x + ce0 = 0
    CI^T x + ci0 >= 0
	 
 The matrix and vectors dimensions are as follows:
     G: n * n
		g0: n
				
		CE: n * p
	 ce0: p
				
	  CI: n * m
   ci0: m

     x: n
 
The function returns true if there is a solution.
If there is, it writes it into x.
NOTE: If x is of the proper size, it is taken as a guess at the solution. 
 
 References: D. Goldfarb, A. Idnani. A numerically stable dual method for solving
             strictly convex quadratic programs. Mathematical Programming 27 (1983) pp. 1-33.

 Notes:
  1. pay attention in setting up the vectors ce0 and ci0. 
	   If the constraints of your problem are specified in the form 
	   A^T x = b and C^T x >= d, then you should set ce0 = -b and ci0 = -d.  
	   2. This takes B (BG = I) and θ (= B g0)
    
 Author: Luca Di Gaspero
  			 DIEGM - University of Udine, Italy
				 l.digaspero@uniud.it
				 http://www.diegm.uniud.it/digaspero/
 
 The author will be grateful if the researchers using this software will
 acknowledge the contribution of this function in their research papers.

LICENSE

This file is part of QuadProg++: a C++ library implementing
the algorithm of Goldfarb and Idnani for the solution of a (convex) 
Quadratic Programming problem by means of an active-set dual method.
Copyright (C) 2007-2009 Luca Di Gaspero.
Copyright (C) 2009 Eric Moyer.  
Copyright 2017 Eyal Minsky-Fenick (much simplified).

QuadProg++ is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

QuadProg++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with QuadProg++. If not, see <http://www.gnu.org/licenses/>.

*/


#include <Eigen/Core>
#include <utility>
namespace ComputationalLib {

  using namespace Eigen;
  extern bool solve_quadprog(const MatrixXd &B, const VectorXd &theta, 
		      const MatrixXd &CE, const VectorXd &ce0,  
		      const MatrixXd &CI, const VectorXd &ci0, 
		      VectorXd &x);
  extern std::pair<MatrixXd, VectorXd> makeBoxConstraints(VectorXd mins,
						   VectorXd maxes);
  extern double qp_solve(const MatrixXd &G, const VectorXd &g0,
			 const MatrixXd &CE, const VectorXd &ce0,
			 const MatrixXd &CI, const VectorXd &ci0,
			 VectorXd &x, double c);
}
