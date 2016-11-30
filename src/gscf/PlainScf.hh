#pragma once
#include "ScfBase.hh"
#include "ScfStateBase.hh"
#include <krims/ParameterMap.hh>

namespace gscf {

/** Default scf state for a plain scf
 *
 * \tparam ProblemMatrix The type of the Problem matrix object.
 * */
template <typename ProblemMatrix, typename OverlapMatrix>
struct PlainScfState : public ScfStateBase<ProblemMatrix, OverlapMatrix, ProblemMatrix> {
  typedef ScfStateBase<ProblemMatrix, OverlapMatrix, ProblemMatrix> base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::diagmat_type diagmat_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  PlainScfState(probmat_type probmat, const overlap_type& overlap_mat)
        : base_type{std::move(probmat), overlap_mat} {}
};

/** \name  Basic and plain SCF algorithm
 *
 * This just keeps running until maxiter is reached. If you want to add
 * a convergence criterion override the is_converged(state) function.
 *
 * ## Control parameters and their default values
 *   - max_iter: Maximum number of iterations. (Default: 100)
 *   - n_eigenpairs: The number of eigenpairs to seek
 *          (Default: linalgwrap::Constanst<size_type>::all)
 *   - eigensolver_params: krims::Parameter map containing the parameters
 *          for the eigensolver. See linalgwrap/eigensystem.hh or the
 *          documentation of your chosen eigensolver for possible values.
 *          Note that this map can be used to *select* an eigensolver as
 *          well.
 *   - max_error_norm:  If the Frobenius norm of the most recent error
 *                      vector/matrix computed by ``calculate_error``
 *                      is below this value, we consider the iteration
 *                      converged (default: 5e-7)
 *
 * \tparam ScfState  The precise scf state type which is available
 *                   to is_converged and all handler functions.
 *                   This also fixes the type of the problem matrix
 *                   and of the overlap matrix. Should be a the class
 *                   PulayDiisScfState or a subclass of it.
 */
template <typename ScfState>
class PlainScf : public ScfBase<ScfState> {
 public:
  typedef ScfBase<ScfState> base_type;
  typedef typename base_type::state_type state_type;

  typedef typename state_type::probmat_type probmat_type;
  typedef typename state_type::overlap_type overlap_type;
  typedef typename state_type::real_type real_type;
  typedef typename state_type::vector_type vector_type;
  // typedef typename state_type::diagmat_type diagmat_type;
  // typedef typename state_type::scalar_type scalar_type;
  typedef typename state_type::matrix_type matrix_type;

  static_assert(
        std::is_base_of<PlainScfState<probmat_type, overlap_type>, ScfState>::value,
        "ScfState needs to be derived off PulayDiisScfState");

  /** \name Constructor */
  //@{
  /** Construct a plain SCF solver with default parameters */
  PlainScf() {}

  /** Construct a plain SCF solver setting the parameters from the map */
  PlainScf(const krims::ParameterMap& map) : PlainScf() { update_control_params(map); }
  //@}

  /** \name Iteration control */
  ///@{
  /** Update control parameters from Parameter map */
  void update_control_params(const krims::ParameterMap& map) {
    base_type::update_control_params(map);
  }

  /** Get the current settings of all internal control parameters and
   *  update the ParameterMap accordingly.
   */
  void get_control_params(krims::ParameterMap& map) const {
    base_type::get_control_params(map);
  }
  ///@}

  /** Implementation of the SolverBase method */
  void solve_state(state_type& state) const override;
};

//
// ---------------------------------------------------
//

template <typename ScfState>
void PlainScf<ScfState>::solve_state(state_type& state) const {
  assert_dbg(!state.is_failed(), krims::ExcInvalidState("Cannot solve a failed state"));

  while (!base_type::convergence_reached(state)) {
    base_type::start_iteration_step(state);

    // Diagonalise the problem matrix
    state.diagonalised_matrix_ptr = state.problem_matrix_ptr;
    base_type::update_eigenpairs(state);
    state.n_mtx_applies() += state.eigenproblem_stats().n_mtx_applies();

    // We free the extra pointer to the problem matrix object here
    // such that the problem_matrix_ptr is the only guy pointing
    // to it. This has the advantage that the update_problem_matrix
    // does not need to do a copy of the object referred to by the
    // problem_matrix_ptr, but it can just update the current one
    // in-place.
    state.diagonalised_matrix_ptr.reset();
    base_type::update_problem_matrix(state);

    // Calculate the new error.
    base_type::update_last_error_norm(state);

    base_type::end_iteration_step(state);
  }
}

}  // gscf
