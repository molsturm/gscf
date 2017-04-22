//
// Copyright (C) 2017 by the gscf authors
//
// This file is part of gscf.
//
// gscf is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// gscf is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with gscf. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include "ScfBaseKeys.hh"
#include "ScfStateBase.hh"
#include <linalgwrap/Base/Solvers.hh>
#include <linalgwrap/EigensystemSolver.hh>
#include <linalgwrap/eigensystem.hh>
#include <linalgwrap/rescue.hh>

namespace gscf {

DefSolverException1(ExcInnerEigensolverFailed, std::string, details,
                    << "The SCF procedure failed to converge, since an inner "
                       "eigensolver failed: "
                    << details);

/** Class to provide building block functions and a common interface
 * for SCF procedures
 *
 * The error in the base implementation is dummy, i.e. it is tuned such
 * that the iteration will run forever (i.e. until max_iterations is reached).
 *
 * To change this behaviour override ``is_converged`` and/or ``calculate_error``
 *
 * \see gscf::ScfStateBase
 * */
template <typename State>
class ScfBase : public linalgwrap::IterativeWrapper<linalgwrap::SolverBase<State>> {
  static_assert(IsScfState<State>::value,
                "State needs to be a type derived from ScfStateBase");

 public:
  typedef linalgwrap::IterativeWrapper<linalgwrap::SolverBase<State>> base_type;
  typedef typename base_type::state_type state_type;

  /** \name Types forwarded from ScfState */
  ///@{
  /** The type of the problem matrix as determined by the traits */
  typedef typename state_type::probmat_type probmat_type;

  /** The type of the overlap matrix as determined by the traits */
  typedef typename state_type::overlap_type overlap_type;

  /** The type of the problem matrix as determined by the traits */
  typedef typename state_type::diagmat_type diagmat_type;

  /** The type of the scalars as determined by the traits */
  typedef typename state_type::scalar_type scalar_type;

  /** The type of the scalars as determined by the traits */
  typedef typename state_type::real_type real_type;

  /** The type of the size indices as determined by the traits */
  typedef typename state_type::size_type size_type;

  /** The type of the stored matrices as determined by the traits */
  typedef typename state_type::matrix_type matrix_type;

  /** The type of the vector as determined by the traits */
  typedef typename state_type::vector_type vector_type;
  ///@}

  /** \name Iteration control */
  ///@{
  /** The number of eigenpairs to calculate in the SCF procedure */
  size_type n_eigenpairs = linalgwrap::Constants<size_type>::all;

  /** The parameters for the inner eigensolver */
  krims::GenMap eigensolver_params;

  /** Maximum value the Frobenius norm of the
   *  most recent error vector/matrix may have. */
  real_type max_error_norm = 5e-7;

  /** Bulk-update control parameters from a parameter map.
   *
   * For the list of available keys, see ScfBaseKeys.hh
   */
  void update_control_params(const krims::GenMap& map) {
    base_type::update_control_params(map);
    n_eigenpairs = map.at(ScfBaseKeys::n_eigenpairs, n_eigenpairs);
    max_error_norm = map.at(ScfBaseKeys::max_error_norm, max_error_norm);
    eigensolver_params = map.submap(ScfBaseKeys::eigensolver_params);
    user_eigensolver_tolerance = map.exists(linalgwrap::EigensystemSolverKeys::tolerance);
  }

  /** Has the user provided tolerance settings for the eigensolver? */
  bool user_eigensolver_tolerance = false;

  /** Get the current settings of all internal control parameters and
   *  update the GenMap accordingly.
   */
  void get_control_params(krims::GenMap& map) const {
    base_type::get_control_params(map);
    map.update(ScfBaseKeys::n_eigenpairs, n_eigenpairs);
    map.update(ScfBaseKeys::max_error_norm, max_error_norm);
    map.update(ScfBaseKeys::eigensolver_params, eigensolver_params);
  }

  /** Is the iteration error already smaller than the
   * required norm for convergence ? */
  bool is_converged(const state_type& s) const override {
    // Norm of the current error is below threshold
    return s.last_error_norm < max_error_norm;
  }
  ///@}

  /** \name Run an SCF
   */
  ///@{
  /** \brief Run the scf solver with the provided problem matrix
   *  and overlap matrix and return the final state.
   *
   *  If the solver does not manage to achieve convergence a
   *  SolverException is thrown as the state's fail bit will be set
   *  accompanied with an appropriate fail message.
   * as well
   **/
  virtual state_type solve(probmat_type probmat_bb,
                           const overlap_type& overlap_bb) const {
    state_type state{std::move(probmat_bb), overlap_bb};
    this->solve_state(state);
    return state;
  }

  /** \brief Run the solver on a problem, starting from a guess state
   *
   * Here we use the ScfStateBase in order to be able to use
   * states of potentially different state_type as well.
   */
  template <typename GuessState,
            typename = krims::enable_if_t<
                  std::is_base_of<ScfStateBase<probmat_type, overlap_type, diagmat_type>,
                                  krims::remove_reference_t<GuessState>>::value>>
  state_type solve_with_guess(probmat_type probmat_bb, const overlap_type& overlap_bb,
                              GuessState&& guess_state) const {
    // Create a new state and install the guess state:
    state_type state{std::move(probmat_bb), overlap_bb};
    state.obtain_guess_from(std::forward<GuessState>(guess_state));
    this->solve_state(state);
    return state;
  }

  /** \brief Run the solver on a problem, starting from a guess state
   *
   * Here we use the EigensolverStateBase in order to be able to use
   * states of potentially different state_type as well.
   */
  template <typename GuessState,
            typename = krims::enable_if_t<
                  std::is_base_of<ScfStateBase<probmat_type, overlap_type, diagmat_type>,
                                  krims::remove_reference_t<GuessState>>::value>>
  state_type solve_with_guess(GuessState&& guess_state) const {
    // Create a new state and install the guess state:
    state_type state{guess_state.problem_matrix(), guess_state.overlap_matrix()};
    state.obtain_guess_from(std::forward<GuessState>(guess_state));
    this->solve_state(state);
    return state;
  }
  ///@}

 protected:
  /** \name Handler functions
   * Various virtual handler functions, which are called when
   * certain events happen.
   *
   * Further handler functions can be found in IterativeSolver and
   * SolverBase of linalgwrap.
   */
  ///@{
  /** Handler which is called once the problem matrix has been diagonalised
   *  and new eigenpairs are obtained for this problem matrix
   */
  virtual void on_update_eigenpairs(state_type&) const {}

  /** Handler which is called once the new problem matrix has been formed
   *  from the current set of eigenvectors
   */
  virtual void on_update_problem_matrix(state_type&) const {}
  ///@}

  /** \name Scf building blocks.
   * Various virtual handler functions, which are called when
   * certain events happen.
   */
  ///@{
  /** Calculate the error from a provided scf state, which should be
   *  a PulayDiisScfState.
   *
   *  This implementation calculates the residual error, i.e. the difference
   *  vectors between the current and the most recent eigenvectors.
   *
   * \note This function is called once the updated problem matrix has
   * been obtained.
   **/
  virtual matrix_type calculate_error(const state_type&) const;

  /** Update the last_error_norm using the calculate_error function. */
  void update_last_error_norm(state_type& s) const {
    s.last_error_norm = norm_frobenius(calculate_error(s));
  }

  /** \brief Solve the eigensystem currently represented by the
   *  diagonalised_matrix_ptr and the overlap_matrix of the SCF state
   *  and place results back inside the state.
   *
   *  After everything is done, we call the on_update_eigenpairs.
   *
   *  Notice, that we just re-assign the shared pointers
   *  and that their magic is in charge for cleaning up
   *  unused storage automatically.
   */
  void update_eigenpairs(state_type& s) const;

  /** \brief Update the problem matrix and replace the current
   *  one in the SCF state.
   *
   *  Notice, that we just re-assign the shared pointers
   *  and that their magic is in charge for cleaning up
   *  unused storage automatically.
   */
  void update_problem_matrix(state_type& s) const;

  /** \brief Compute a sensible tolerance value to use for
   * inner (eigen)solvers in this iteration */
  real_type inner_solver_tolerance(state_type& s) const {
    const real_type fac = 0.01;

    // TODO really do some playing around here!
    return std::numeric_limits<real_type>::epsilon();

    // For the first iteration there is no last error
    real_type tolerance =
          (s.n_iter() == 1) ? fac * max_error_norm : fac * s.last_error_norm;
    assert_finite(tolerance);

    // No undershooting of numeric epsilon
    return std::max(tolerance, std::numeric_limits<real_type>::epsilon());
  }
  ///@}
};

//
// ------------------------------------
//

template <typename ScfState>
typename ScfBase<ScfState>::matrix_type ScfBase<ScfState>::calculate_error(
      const state_type& s) const {
  if (s.previous_eigensolution().evectors().n_vectors() == 0) {
    // Cannot compute any matrix, since no previous vectors
    // to compare against:
    return matrix_type{{linalgwrap::Constants<scalar_type>::invalid}};
  }

  const auto& prev_evec = s.previous_eigensolution().evectors();
  const auto& cur_evec = s.eigensolution().evectors();
  matrix_type ret(prev_evec.n_elem(), prev_evec.n_vectors(), false);
  for (size_type j = 0; j < cur_evec.n_vectors(); ++j) {
    for (size_type i = 0; i < cur_evec.n_elem(); ++i) {
      const vector_type& c = cur_evec[j];
      const vector_type& p = prev_evec[j];
      ret(i, j) = c(i) - p(i);
    }  // i == row
  }    // j == col

  return ret;
}

template <typename ScfState>
void ScfBase<ScfState>::update_eigenpairs(state_type& s) const {
  using namespace linalgwrap;

  assert_dbg(s.diagonalised_matrix_ptr != nullptr,
             krims::ExcInvalidState("update_eigenpairs needs a valid "
                                    "matrix pointer inside "
                                    "s.diagonalised_matrix_ptr"));

  // TODO make use of SCF tolerance somehow

  // Matrix and eigenproblem setup:
  typedef Eigenproblem</* herm= */ true, diagmat_type, overlap_type> eprob_type;
  eprob_type problem(s.diagonalised_matrix(), s.overlap_matrix(), n_eigenpairs);

  // Setup solver state with the eigensolution from the previous run.
  EigensystemSolverState<eprob_type> state(problem);
  state.obtain_guess_from(s.eigensolution());

  // Setup the solver. If the user did not supply a tolerance,
  // then we do.
  EigensystemSolver<eprob_type> solver{eigensolver_params};
  if (!user_eigensolver_tolerance) {
    solver.tolerance = inner_solver_tolerance(s);
  }

  try {
    solver.solve_state(state);
    s.push_new_eigensolution(std::move(state.eigensolution()),
                             {state.n_iter(), state.n_mtx_applies()});
  } catch (const linalgwrap::SolverException& e) {
    linalgwrap::rescue::failed_eigenproblem(problem, eigensolver_params);
    solver_assert(false, s, ExcInnerEigensolverFailed(e.extra()));
  }
  // Call the handler
  on_update_eigenpairs(s);
}

template <typename ScfState>
void ScfBase<ScfState>::update_problem_matrix(state_type& s) const {
  if (!s.problem_matrix_ptr.unique()) {
    // s is not the only thing referencing the object behind the
    // problem_matrix_ptr shared pointer, so we need to copy the object
    // first. This assures that objects that depend on the SCF state (i.e.
    // those which actually copyied the shared pointer to some internal
    // place) still access the old, non-updated problem matrix.
    //
    // A typical scenario for this could be algorithms similar to the
    // Pulay-DIIS, where a linear combination Problem matrices are used
    // to produce the next guess for the eigenproblem.
    // There the integrity of the problem matrices needs to be preserved,
    // as we need them in the future.

    // Copy the old problem matrix:
    auto new_problem_matrix_ptr = std::make_shared<probmat_type>(*s.problem_matrix_ptr);

    // Replace the current one in the state:
    s.problem_matrix_ptr = std::move(new_problem_matrix_ptr);
  }

  // Obtain the expected update key from the problem matrix
  // and update the problem matrix:
  const std::string key = s.problem_matrix().scf_update_key();
  const auto const_evec_ptr = s.eigensolution().evectors_ptr;
  s.problem_matrix().update({{key, const_evec_ptr}});

  on_update_problem_matrix(s);
}

}  // namespace gscf
