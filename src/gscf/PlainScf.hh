#pragma once
#include "ScfBase.hh"
#include "ScfStateBase.hh"
#include <krims/ParameterMap.hh>

namespace gscf {

/** Default scf state for a plain scf
 *
 * \tparam ProblemMatrix The type of the Problem matrix object.
 * */
template <typename ProblemMatrix>
struct PlainScfState : public ScfStateBase<ProblemMatrix, ProblemMatrix> {
  typedef ScfStateBase<ProblemMatrix, ProblemMatrix> base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::diagmat_type diagmat_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  PlainScfState(probmat_type probmat, const matrix_type& overlap_mat)
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
 *
 * \tparam ProblemMatrix  The type of the problem matrix to be solved
 * \tparam ScfState       The precise scf state type which is available
 *                        to is_converged and all handler functions.
 */
template <typename ProblemMatrix,
          typename ScfState = PlainScfState<ProblemMatrix>>
class PlainScf : public ScfBase<ScfState> {
public:
  typedef ScfBase<ScfState> base_type;
  typedef typename base_type::state_type state_type;

  typedef typename state_type::probmat_type probmat_type;
  // typedef typename state_type::diagmat_type diagmat_type;
  // typedef typename state_type::scalar_type scalar_type;
  // typedef typename state_type::size_type size_type;
  typedef typename state_type::matrix_type matrix_type;

  static_assert(std::is_same<ProblemMatrix, probmat_type>::value,
                "The ProblemMatrix type specified and the one implicit in the "
                "ScfState class have to agree.");

  /** \name Constructor */
  //@{
  /** Construct a plain SCF solver with default parameters */
  PlainScf() {}

  /** Construct a plain SCF solver setting the parameters from the map */
  PlainScf(const krims::ParameterMap& map) : PlainScf() {
    update_control_params(map);
  }
  //@}

  /** \name Iteration control */
  ///@{
  /** Update control parameters from Parameter map */
  void update_control_params(const krims::ParameterMap& map) {
    base_type::update_control_params(map);
  }
  ///@}

  /** Implementation of the SolverBase method */
  void solve_state(state_type& state) const override;
};

template <typename ProblemMatrix, typename ScfState>
void PlainScf<ProblemMatrix, ScfState>::solve_state(state_type& state) const {
  while (!base_type::convergence_reached(state)) {
    base_type::start_iteration_step(state);

    // Diagonalise the problem matrix
    state.diagonalised_matrix_ptr() = state.problem_matrix_ptr();
    base_type::update_eigenpairs(state);

    // We free the extra pointer to the problem matrix object here
    // such that the problem_matrix_ptr is the only guy pointing
    // to it. This has the advantage that the update_problem_matrix
    // does not need to do a copy of the object referred to by the
    // problem_matrix_ptr, but it can just update the current one
    // in-place.
    state.diagonalised_matrix_ptr().reset();
    base_type::update_problem_matrix(state);

    base_type::end_iteration_step(state);
  }
}

}  // gscf
