/* 
 
 Author: Luca Di Gaspero
 DIEGM - University of Udine, Italy
 l.digaspero@uniud.it
 http://www.diegm.uniud.it/digaspero/
 
 Rewritten by Eyal Minsky-Fenick

 LICENSE 
 
 This file is part of QuadProg++: a C++ library implementing
 the algorithm of Goldfarb and Idnani for the solution of a (convex) 
 Quadratic Programming problem by means of an active-set dual method.
 Copyright (C) 2007-2009 Luca Di Gaspero.
 Copyright (C) 2009 Eric Moyer.  
 
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

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "qp.hpp"
#include <vector>

#ifndef EPSILON
#define EPSILON .0000000000000002220446049250313080847263336181640625
#endif

using namespace Eigen;
using std::vector;

namespace linreg {

  // Utility functions for updating some data needed by the solution method
  inline void update_r(const MatrixXd &R, VectorXd &r,
		       const VectorXd &d, int iq) {
    double sum;
    
    /* setting of r = R^-1 d */
    for(int i = iq - 1; i >= 0; i--) {
      sum = 0.0;
      for(int j = i + 1; j < iq; j++)
	sum += R(i, j) * r(j);
      r(i) = (d(i) - sum) / R(i, i);
    }
  }
  
  bool add_constraint(MatrixXd &R, MatrixXd &J, VectorXd &d,
		      unsigned &iq, double &rnorm);
  void delete_constraint(MatrixXd &R, MatrixXd &J, VectorXi &A, VectorXd &u,
			 unsigned n, unsigned p, unsigned &iq, unsigned l);
  
  // Utility functions for computing the scalar product and the euclidean 
  // distance between two numbers
  inline double distance(double a, double b) {
    double a1, b1, t;
    a1 = fabs(a);
    b1 = fabs(b);
    if (a1 > b1) {
      t = (b1 / a1);
      return a1 * sqrt(1.0 + t * t);
    } else
      if (b1 > a1) {
	t = (a1 / b1);
	return b1 * sqrt(1.0 + t * t);
      }
    return a1 * sqrt(2.0);
  } 
  
  
  // The Solving function, implementing the Goldfarb-Idnani method
  
  bool solve_quadprog(const MatrixXd &B, const VectorXd &theta, 
		      const MatrixXd &CE, const VectorXd &ce0,  
		      const MatrixXd &CI, const VectorXd &ci0, 
		      VectorXd &x) {
    unsigned n = B.cols(), p = ce0.size(), m = ci0.size(), l;
    if(B.rows() != B.cols()) {
      std::ostringstream msg;
      msg << "The matrix B is not square (" << B.rows()
	  << 'x' << B.cols() << ')';
      throw std::logic_error(msg.str());
    }
    if(p) {
      if(CE.rows() != B.rows()) {
	std::ostringstream msg;
	msg << "The matrix CE is incompatible (number of rows: " 
	    << CE.rows() << "; expecting: " << B.rows() << ')';
	throw std::logic_error(msg.str());
      }
      if(p != CE.cols()) {
	std::ostringstream msg;
	msg << "The vector ce0 is incompatible (size: " 
	    << ce0.size() << "; expecting: " << CE.cols() << ')';
	throw std::logic_error(msg.str());
      }
    }
    if(m) {
      if(CI.rows() != B.rows()) {
	std::ostringstream msg;
	msg << "The matrix CI is incompatible (number of rows: " 
	    << CI.rows() << "; expecting: " << B.rows() << ')';
	throw std::logic_error(msg.str());
      }
      if(m != CI.cols()) {
	std::ostringstream msg;
	msg << "The vector ci0 is incompatible (size: " 
	    << ci0.size() << "; expecting: " << CI.cols() << ')';
	throw std::logic_error(msg.str());
      }
    }
    
    if(x.size() != n)
      x = theta;
    
    int ip; // this is the index of the constraint to be added to the active set
    MatrixXd R = MatrixXd::Zero(n, n), J = B;
    VectorXd s(m + p), z(n), r(m + p), d(n), np(n), u(m + p), x_old(n), u_old(m + p);
    double psi, c1, c2, ss, R_norm;
    double t, t1, t2; /* t is the step lenght, which is the minimum of the partial step length t1 and the full step length t2 */
    VectorXi A(m + p), A_old(m + p), iai(m + p);
    unsigned q, iq;
    vector<bool> iaexcl(m + p);
    
    /* p is the number of equality constraints */
    /* m is the number of inequality constraints */
    q = 0;  /* size of the active set A (containing the indices of the active constraints) */     
    
    /* Add equality constraints to the working set A */
    iq = 0;
    for(unsigned i = 0; i < p; i++) {
      np = CE.col(i);
      d = J * np;
      unsigned j = J.cols() - iq;
      z = J.rightCols(j) * d.tail(j);
      update_r(R, r, d, iq);
      
      /* compute full step length t2: i.e., the minimum step in primal space s.t. the contraint 
	 becomes feasible */
      t2 = 0.0;
      if(z.dot(z) > EPSILON) // i.e. z != 0
	t2 = (-z.dot(x) - ce0(i)) / z.dot(np);
      
      /* set x = x + t2 * z */
      x += z * t2;
      
      /* set u = u+ */
      u(iq) = t2;
      u.head(iq) -= r.head(iq) * t2;
      
      /* compute the new solution value */
      A(i) = -i - 1;
       
      if(!add_constraint(R, J, d, iq, R_norm))
	// Equality constraints are linearly dependent
	throw std::runtime_error("Constraints are linearly dependent");
    }
    
    /* set iai = K \ A */
    for(unsigned i = 0; i < m; i++)
      iai(i) = i;
    
    while(true) {
      /* step 1: choose a violated constraint */
      for (unsigned i = p; i < iq; i++)
	iai(A(i)) = -1;
      
      /* compute s(x) = ci^T * x + ci0 for all elements of K \ A */
      ss = 0.0;
      psi = 0.0; /* this value will contain the sum of all infeasibilities */
      ip = 0; /* ip will be the index of the chosen violated constraint */
      for(unsigned i = 0; i < m; i++) {
	iaexcl[i] = true;
	s(i) = CI.col(i).dot(x) + ci0(i);
	psi += std::min(0.0, s(i));
      }
      if (fabs(psi) <= m * EPSILON * c1 * c2 * 100.0)
	/* numerically there are not infeasibilities anymore */
	return(true);
      
      /* save old values for u and A */
      u_old = u;
      A_old = A;
      /* and for x */
      x_old = x;
      
    l2: /* Step 2: check for feasibility and determine a new S-pair */
      for(unsigned i = 0; i < m; i++) {
	if(s(i) < ss && ~iai(i) && iaexcl[i]) {
	  ss = s(i);
	  ip = i;
	}
      }
      if(ss >= 0)
	return(true);
      
      /* set np = n(ip) */
      np = CI.col(ip);
      /* set u = (u 0)^T */
      u(iq) = 0;
      /* add ip to the active set A */
      A(iq) = ip;
      
    l2a:/* Step 2a: determine step direction */
      /* compute z = H np: the step direction in the primal space (through J, see the paper) */
      d = J * np;
      z = J.rightCols(J.cols() - iq) * d.tail(J.cols() - iq);
      /* compute N* np (if q > 0): the negative of the step direction in the dual space */
      update_r(R, r, d, iq);
      
      /* Step 2b: compute step length */
      l = 0;
      /* Compute t1: partial step length (maximum step in dual space without violating dual feasibility */
      t1 = INFINITY; /* +inf */
      /* find the index l s.t. it reaches the minimum of u+(x) / r */
      for(unsigned k = p; k < iq; k++) {
	if(r(k) > 0 && u(k) / r(k) < t1) {
	  t1 = u(k) / r(k);
	  l = A(k);
	}
      }
      /* Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible */
      if(z.dot(z) > EPSILON) // i.e. z != 0
	t2 = -s(ip) / z.dot(np);
      else
	t2 = INFINITY; /* +inf */
      
      /* the step is chosen as the minimum of t1 and t2 */
      t = std::min(t1, t2);
      
      /* Step 2c: determine new S-pair and take step: */
      
      /* case (i): no step in primal or dual space */
      if (t == INFINITY)
	return(false);
      /* case (ii): step in dual space */
      if (t2 == INFINITY) {
	// set u = u +  t * (-r 1) and drop constraint l from the active set A
	u.head(iq) -= r.head(iq) * t;
	u(iq) += t;
	iai(l) = l;
	delete_constraint(R, J, A, u, n, p, iq, l);
	goto l2a;
      }
      
      /* case (iii): step in primal and dual space */
      
      x += z * t;
      /* u = u + t * (-r 1) */
      u.head(iq) -= r.head(iq) * t;
      u(iq) += t;
      
      if(fabs(t - t2) < EPSILON) {
	/* full step has taken */
	/* add constraint ip to the active set*/
	if(!add_constraint(R, J, d, iq, R_norm)) {
	  iaexcl[ip] = false;
	  delete_constraint(R, J, A, u, n, p, iq, ip);
	  for(unsigned i = 0; i < m; i++)
	    iai(i) = i;
	  for(unsigned i = p; i < iq; i++) {
	    A(i) = A_old(i);
	    u(i) = u_old(i);
	    iai(A(i)) = -1;
	  }
	  x = x_old;
	  goto l2; /* go to step 2 */
	}
	iai(ip) = -1;
	continue;
      }
      
      /* a patial step has taken */
      /* drop constraint l */
      iai(l) = l;
      delete_constraint(R, J, A, u, n, p, iq, l);
      
      /* update s(ip) = CI * x + ci0 */
      s(ip) = CI.col(ip).dot(x) + ci0(ip);
      
      goto l2a;
    }
  }
   
  bool add_constraint(MatrixXd &R, MatrixXd &J, VectorXd &d,
		      unsigned &iq, double &R_norm) {
    unsigned n = d.size();
    double cc, ss, h, t1, t2, xny;
    
    /* we have to find the Givens rotation which will reduce the element
       d(j) to zero.
       if it is already zero we don't have to do anything, except of
       decreasing j */  
    for(unsigned j = n - 1; j >= iq + 1; j--) {
      /* The Givens rotation is done with the ublas::matrix (cc cs, cs -cc).
	 If cc is one, then element (j) of d is zero compared with element
	 (j - 1). Hence we don't have to do anything. 
	 If cc is zero, then we just have to switch columns (j) and (j - 1) 
	 of J. Since we only switch columns in J, we have to be careful how we
	 update d depending on the sign of gs.
	 Otherwise we have to apply the Givens rotation to these columns.
	 The i - 1 element of d has to be updated to h. */
      cc = d(j - 1);
      ss = d(j);
      h = distance(cc, ss);
      if (fabs(h) < EPSILON) // h == 0
	continue;
      d(j) = 0.0;
      ss = ss / h;
      cc = cc / h;
      if (cc < 0.0) {
	cc = -cc;
	ss = -ss;
	d(j - 1) = -h;
      } else
       d(j - 1) = h;
      xny = ss / (1.0 + cc);
      for(unsigned k = 0; k < n; k++) {
	t1 = J(k, j - 1);
	t2 = J(k, j);
	J(k, j - 1) = t1 * cc + t2 * ss;
	J(k, j) = xny * (t1 + J(k, j - 1)) - t2;
      }
    }
    /* update the number of constraints added*/
    iq++;
    /* To update R we have to put the iq components of the d vector
      into column iq - 1 of R
    */
    
    R.col(iq - 1).head(iq) = d.head(iq);
    
    if (fabs(d(iq - 1)) <= EPSILON * R_norm) {
      // problem degenerate
      return false;
    }
    R_norm = std::max<double>(R_norm, fabs(d(iq - 1)));
    return true;
 }
  
  void delete_constraint(MatrixXd &R, MatrixXd &J, VectorXi &A, VectorXd &u,
			 unsigned n, unsigned p, unsigned &iq, unsigned l) {
    unsigned qq;
    double cc, ss, h, xny, t1, t2;
    
   /* Find the index qq for active constraint l to be removed */
    for(qq = p; qq < iq && A(qq) != l; qq++);
    
    /* remove the constraint from the active set and the duals */
    for (unsigned i = qq; i < iq - 1; i++) {
      A(i) = A(i + 1);
      u(i) = u(i + 1);
      R.col(i) = R.col(i + 1);
    }
    
    A(iq - 1) = A(iq);
    u(iq - 1) = u(iq);
    A(iq) = 0;
    u(iq) = 0;
    R.col(iq - 1).head(iq).setZero();
    /* constraint has been fully removed */
    if(!--iq)
      return;
    
    for(unsigned j = qq; j < iq; j++) {
      cc = R(j, j);
      ss = R(j + 1, j);
      h = distance(cc, ss);
      if (fabs(h) < EPSILON) // h == 0
	continue;
      cc = cc / h;
      ss = ss / h;
      R(j + 1, j) = 0;
      if(cc < 0) {
	R(j, j) = -h;
	cc = -cc;
	ss = -ss;
      }
      else
	R(j, j) = h;
      
      xny = ss / (1.0 + cc);
      for(unsigned k = j + 1; k < iq; k++) {
	t1 = R(j, k);
	t2 = R(j + 1, k);
	R(j, k) = t1 * cc + t2 * ss;
	R(j + 1, k) = xny * (t1 + R(j, k)) - t2;
	t1 = J(k, j);
	t2 = J(k, j + 1);
	J(k, j) = t1 * cc + t2 * ss;
	J(k, j + 1) = xny * (J(k, j) + t1) - t2;
      }
    }
  }

  std::pair<MatrixXd, VectorXd> makeBoxConstraints(VectorXd &mins,
						   VectorXd &maxes) {
    int i;
    for(int j = 0; j < mins.size(); j++) {
      if(maxes(j) < INFINITY)
	i++;
      if(mins(j) > -INFINITY)
	i++;
    }

    MatrixXd ci(mins.size(), i);
    VectorXd ci0(i);
    
    for(int j = 0; j < mins.size(); j++) {
      if(mins(j) > -INFINITY) {
	ci.col(--i).setZero();
	ci(j, i) = 1;
	ci0(i) = -mins(j);
      }
      if(maxes(j) < INFINITY) {
	ci.col(--i).setZero();
	ci(j, i) = -1;
	ci0(i) = maxes(j);
      }
    }
    return(std::make_pair(ci, ci0));
  }
    
  
}
