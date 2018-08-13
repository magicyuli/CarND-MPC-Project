#include "MPC.h"

#include <iostream>

#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>

#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

const size_t N = 30;
const double dt = 0.05;
const size_t latency_offset = 0.1 / dt;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

const double ref_v = 40.0;

const size_t x_start = 0;
const size_t y_start = x_start + N;
const size_t v_start = y_start + N;
const size_t psi_start = v_start + N;
const size_t cte_start = psi_start + N;
const size_t e_psi_start = cte_start + N;
const size_t a_start = e_psi_start + N;
const size_t delta_start = a_start + N - 1;

// 6 variables and 2 actuators
// x, y, psi, v, cte, e_psi
const size_t n_vars = 6 * N + 2 * (N - 1);
const size_t n_constraints = 6 * N;

class FG_eval {
public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;

  FG_eval(Eigen::VectorXd coeffs) {
    this->coeffs = coeffs;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

  void operator()(ADvector &fg, const ADvector &vars) {
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    // calculate error f(x)
    fg[0] = 0;
    for (size_t t = 0; t < N; t++) {
      fg[0] += CppAD::pow(vars[cte_start + t], 2);
      fg[0] += CppAD::pow(vars[e_psi_start + t], 2);
      fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);
    }
    for (size_t t = 0; t < N - 1; t++) {
      fg[0] += CppAD::pow(vars[a_start + t], 2) * 3e2;
      fg[0] += CppAD::pow(vars[delta_start + t], 2) * 3e2;
    }
    for (size_t t = 0; t < N - 2; t++) {
      fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2) * 1e3;
      fg[0] += CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2) * 1e3;
    }

    // calculate constraints g(x)
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + e_psi_start] = vars[e_psi_start];
    for (size_t t = 0; t < N - 1; t++) {
      // at t
      auto x_t = vars[x_start + t];
      auto y_t = vars[y_start + t];
      auto v_t = vars[v_start + t];
      auto psi_t = vars[psi_start + t];
      auto a_t = vars[a_start + t];
      auto delta_t = vars[delta_start + t];

      // at t + 1
      auto x_tt = vars[x_start + t + 1];
      auto y_tt = vars[y_start + t + 1];
      auto v_tt = vars[v_start + t + 1];
      auto psi_tt = vars[psi_start + t + 1];
      auto cte_tt = vars[cte_start + t + 1];
      auto e_psi_tt = vars[e_psi_start + t + 1];

      fg[1 + x_start + t + 1] = x_tt - (x_t + v_t * CppAD::cos(psi_t) * dt);
      fg[1 + y_start + t + 1] = y_tt - (y_t + v_t * CppAD::sin(psi_t) * dt);
      fg[1 + v_start + t + 1] = v_tt - (v_t + a_t * dt);
      fg[1 + psi_start + t + 1] = psi_tt - (psi_t + v_t / Lf * delta_t * dt);
      fg[1 + cte_start + t + 1] =
        cte_tt -
        (coeffs[0] + coeffs[1] * x_t + coeffs[2] * CppAD::pow(x_t, 2) + coeffs[3] * CppAD::pow(x_t, 3) - y_t +
         v_t * CppAD::sin(psi_t) * dt);
      fg[1 + e_psi_start + t + 1] =
        e_psi_tt -
        (psi_t - CppAD::atan(coeffs[1] + 2 * coeffs[2] * x_t + 3 * coeffs[3] * CppAD::pow(x_t, 2)) +
         v_t / Lf * delta_t * dt);
    }
    std::cout << "fg: " << fg << std::endl;
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}

MPC::~MPC() {}

Solution MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  typedef CPPAD_TESTVECTOR(
  double) Dvector;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (size_t i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }
  vars[x_start] = state[0];
  vars[y_start] = state[1];
  vars[v_start] = state[2];
  vars[psi_start] = state[3];
  vars[cte_start] = state[4];
  vars[e_psi_start] = state[5];

  // Set lower and upper limits for variables.
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  for (size_t i = 0; i < a_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  // simulate latency by fixing a for the latency period
  for (size_t i = a_start; i < a_start + latency_offset; i++) {
    vars_lowerbound[i] = prev_a;
    vars_upperbound[i] = prev_a;
  }
  for (size_t i = a_start + latency_offset; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
  // simulate latency by fixing delta for the latency period
  for (size_t i = delta_start; i < delta_start + latency_offset; i++) {
    vars_lowerbound[i] = prev_delta;
    vars_upperbound[i] = prev_delta;
  }
  for (size_t i = delta_start + latency_offset; i < n_vars; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (size_t i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0.0;
    constraints_upperbound[i] = 0.0;
  }
  constraints_lowerbound[x_start] = state[0];
  constraints_lowerbound[y_start] = state[1];
  constraints_lowerbound[v_start] = state[2];
  constraints_lowerbound[psi_start] = state[3];
  constraints_lowerbound[cte_start] = state[4];
  constraints_lowerbound[e_psi_start] = state[5];

  constraints_upperbound[x_start] = state[0];
  constraints_upperbound[y_start] = state[1];
  constraints_upperbound[v_start] = state[2];
  constraints_upperbound[psi_start] = state[3];
  constraints_upperbound[cte_start] = state[4];
  constraints_upperbound[e_psi_start] = state[5];

  // object that computes objective and constraints
  FG_eval fg_eval(std::move(coeffs));

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result <Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
    options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
    constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  std::cout << "ok: " << (solution.status == CppAD::ipopt::solve_result<Dvector>::success) << std::endl;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost: " << cost << std::endl;

  Solution res{};
  for (size_t t = 0; t < N; t++) {
    res.x.push_back(solution.x[x_start + t]);
    res.y.push_back(solution.x[y_start + t]);
  }
  res.a = solution.x[a_start + latency_offset];
  res.delta = solution.x[delta_start + latency_offset];

  prev_a = res.a;
  prev_delta = res.delta;

  return res;
}
