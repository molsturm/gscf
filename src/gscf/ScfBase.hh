#pragma once
#include "ScfBaseKeys.hh"
#include "ScfStateBase.hh"
#include <linalgwrap/Base/Solvers.hh>
#include <linalgwrap/EigensystemSolver.hh>
#include <linalgwrap/eigensystem.hh>

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
class ScfBase
      : public linalgwrap::IterativeWrapper<linalgwrap::SolverBase<State>> {
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
  ///@}

  /** \name Iteration control */
  ///@{
  /** The number of eigenpairs to calculate in the SCF procedure */
  size_type n_eigenpairs = linalgwrap::Constants<size_type>::all;

  /** The parameters for the inner eigensolver */
  krims::ParameterMap eigensolver_params;

  /** Maximum value the Frobenius norm of the
   *  most recent error vector/matrix may have. */
  real_type max_error_norm = 5e-7;

  /** Bulk-update control parameters from a parameter map.
   *
   * For the list of available keys, see ScfBaseKeys.hh
   */
  void update_control_params(const krims::ParameterMap& map) {
    base_type::update_control_params(map);
    n_eigenpairs = map.at(ScfBaseKeys::n_eigenpairs, n_eigenpairs);
    max_error_norm = map.at(ScfBaseKeys::max_error_norm, max_error_norm);
    eigensolver_params = map.submap(ScfBaseKeys::eigensolver_params);
    user_eigensolver_tolerance =
          map.exists(linalgwrap::EigensystemSolverKeys::tolerance);
  }

  /** Has the user provided tolerance settings for the eigensolver? */
  bool user_eigensolver_tolerance = false;

  /** Get the current settings of all internal control parameters and
   *  update the ParameterMap accordingly.
   */
  void get_control_params(krims::ParameterMap& map) const {
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
   *  SolverException is thnown an the state's fail bit will be set
   *  accompanied with an appropriate fail message.
   * as well
   **/
  virtual state_type solve(probmat_type probmat_bb,
                           const overlap_type& overlapmat_bb) const {
    state_type state{std::move(probmat_bb), overlapmat_bb};
    this->solve_state(state);
    return state;
  }

  /** \brief Run the solver starting from an old SCF state.
   *
   * It is assumed, that the input state is not failed.
   * Note that the fail bit can be unset using the clear_failed() function
   * in order to continue off a failed state using different solver
   * control parameters or methods.
   */
  virtual state_type solve(const state_type& old_state) const {
    assert_dbg(!old_state.is_failed(),
               krims::ExcInvalidState("Cannot make use of a failed state"));
    state_type state{old_state};
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
  /** Calculate the error matrix of a provided scf state.
   *
   *  This implementation calculates the residual error, i.e. the difference
   *  vectors between the current and the most recent eigenvectors.
   **/
  virtual matrix_type calculate_error(const state_type&) const {
    return matrix_type{{linalgwrap::Constants<scalar_type>::invalid}};
  }

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
  ///@}
};

//
// ------------------------------------
//

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
  EigensystemSolverState<eprob_type> state(std::move(problem));
  state.obtain_guess_from(s.eigensolution());

  /*
  if (false && !user_eigensolver_tolerance) {
    // TODO This is absolutely empirical right now.
    //      Check with some literature
    const real_type tolerance = std::max(
          s.last_error_norm / 100., std::numeric_limits<real_type>::epsilon());
    eigensolver_params.update(EigensystemSolverKeys::tolerance, tolerance);
  }
  */

  try {
    // Solve state with parameters:
    EigensystemSolver<eprob_type>{eigensolver_params}.solve_state(state);

    // Update our state:
    s.eigensolution() = state.eigensolution();
    s.n_eigenproblem_iter() = state.n_iter();
  } catch (const linalgwrap::SolverException& e) {
#ifdef DEBUG
    try {
      if (problem.dim() < 1000) {
        std::cerr << "The inner eigensolver failed";

        // Overwrite some user parameters:
        krims::ParameterMap copy(eigensolver_params);
        copy.update(EigensystemSolverKeys::method, "auto");
        copy.update(EigensystemSolverKeys::tolerance,
                    Constants<real_type>::default_tolerance);

        // Solve for full spectrum:
        const auto all = linalgwrap::Constants<size_t>::all;
        auto solresq =
              eigensystem_hermitian(problem.A(), problem.B(), all, copy);

        std::cerr << "  ...  full eigenspectrum of problem:" << std::endl
                  << std::endl;
        std::cerr << "       ";
        std::ostream_iterator<scalar_type> out_it(std::cerr, "  ");
        std::copy(std::begin(solresq.evalues()), std::end(solresq.evalues()),
                  out_it);
        std::cerr << std::endl << std::endl;
      }
    } catch (...) {
      std::cerr << "   ...  but the attempt to solve for the full spectrum "
                   "failed as well."
                << std::endl;
    }
#endif  // DEBUG
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
    auto new_problem_matrix_ptr =
          std::make_shared<probmat_type>(*s.problem_matrix_ptr);

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
