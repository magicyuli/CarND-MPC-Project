#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

struct Solution {
  vector<double> x;
  vector<double> y;
  double a;
  double delta;
};

class MPC {
public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuations.
  Solution Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

private:
  double prev_a = 0.0;
  double prev_delta = 0.0;
};

#endif /* MPC_H */
