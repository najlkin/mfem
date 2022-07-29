#ifndef MATERIALS_SOLVERS_HPP
#define MATERIALS_SOLVERS_HPP

#include "mfem.hpp"

namespace mfem {
namespace materials {

/// Solver for the SPDE method based on a rational approximation with the AAA
/// algorithm. The SPDE method is described in the paper
/// Lindgren, F., Rue, H., Lindström, J. (2011). An explicit link between
/// Gaussian fields and Gaussian Markov random fields: the stochastic partial
/// differential equation approach. Journal of the Royal Statistical Society:
/// Series B (Statistical Methodology), 73(4), 423–498.
/// https://doi.org/10.1111/j.1467-9868.2011.00777.x
///
/// The solver solves the SPDE problem defined as
/// (A)^-\alpha u = b
/// where A is
/// A = div ( Theta(x) grad + Id ) u(x)
/// and \alpha is given as
/// \alpha = (2 nu + dim) / 2.
/// Theta (anisotropy tensor) and nu (smoothness) can be specified in the
/// constructor. Traditionally, the SPDE method requires the specification of
/// a white noise right hands side. SPDESolver accepts arbitrary right hand
/// sides but the solver has only been tested with white noise.
class SPDESolver {
public:
  /// Constructor.
  /// @param diff_coefficient The diffusion coefficient \Theta.
  /// @param nu The coefficient nu, smoothness of the solution.
  /// @param ess_tdof_list Boundary conditions.
  /// @param fespace Finite element space.
  SPDESolver(MatrixConstantCoefficient &diff_coefficient, double nu,
             const Array<int> &ess_tdof_list, ParFiniteElementSpace *fespace);

  /// Solve the SPDE for a given right hand side b. May alter b if
  /// the exponent (alpha) is larger than 1. We avoid copying be default. If you
  /// need b later on, make a copy of it before calling this function.
  void Solve(LinearForm &b, GridFunction &x);

private:
  /// The rational approximation of the SPDE results in multiple
  /// reactio-diffusion PDEs that need to be solved. This call solves the PDE
  /// (div \Theta grad + \alpha I)^exponent x = \beta b.
  void Solve(const LinearForm &b, GridFunction &x, double alpha, double beta,
             int exponent = 1);

  // Each PDE gives rise to a linear system. This call solves the linear system
  // with PCG and Boomer AMG preconditioner.
  void SolveLinearSystem();

  /// Activate repeated solve capabilities. E.g. if the PDE is of the form
  /// A^N x = b. This method solves the PDE A x = b for the first time, and
  /// then uses the solution as RHS for the next solve and so forth.
  void ActivateRepeatedSolve() { repeated_solve_ = true; }

  /// Single solve only.
  void DeactivateRepeatedSolve() { repeated_solve_ = false; }

  /// Writes the solution of the PDE from the previous call to Solve() to the
  /// linear from b (with appropriate transformations).
  void UpdateRHS(LinearForm &b);

  // Compute the coefficients for the rational approximation of the solution.
  void ComputeRationalCoefficients(double exponent);

  // Bilinear forms and corresponding matrices for the solver.
  ParBilinearForm k_;
  ParBilinearForm m_;
  HypreParMatrix stiffness_;
  HypreParMatrix mass_bc_;
  HypreParMatrix mass_0_;

  // Transformation matrices (needed to construct the linear systems and
  // solutions)
  const SparseMatrix *restriction_matrix_;
  const Operator *prolongation_matrix_;

  // Members to solve the linear system.
  Vector X_;
  Vector B_;
  HypreParMatrix *Op_;

  // Information of the finite element space.
  const Array<int> &ess_tdof_list_;
  ParFiniteElementSpace *fespace_ptr_;

  // Coefficients for the rational approximation of the solution.
  Array<double> coeffs_;
  Array<double> poles_;

  // Exponents of the operator
  double nu_ = 0.0;
  double alpha_ = 0.0;
  int integer_order_of_exponent_ = 0;

  // Member to switch to repeated solve capabilities.
  bool repeated_solve_ = false;
  bool integer_order_ = false;
};

} // namespace materials
} // namespace mfem

#endif // MATERIALS_SOLVERS_HPP
