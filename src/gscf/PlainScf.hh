#pragma once
#include "IterationConstants.hh"
#include "ScfBase.hh"
#include "ScfControlBase.hh"
#include "ScfStateBase.hh"
#include <linalgwrap/Constants.hh>
#include <linalgwrap/ParameterMap.hh>

namespace gscf {

/** Default scf state for a plain scf
 *
 * \tparam ProblemMatrix The type of the Problem matrix object.
 * */
template <typename ProblemMatrix>
struct PlainScfState : public ScfStateBase<ProblemMatrix> {
  typedef ScfStateBase<ProblemMatrix> base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  PlainScfState(probmat_type probmat, const matrix_type& overlap_mat)
        : base_type{std::move(probmat), overlap_mat} {}
};

/** \brief Class to provide the basic set of parameters and the convergence
 *  checks the ScfBase class needs.
 *
 *  Consider overriding the is_converged function to provide an actual
 *  convergence check. This version just returns false, i.e. never
 *  converges. The signature is
 *  ```
 *  bool is_converged(const scf_state_type&) const;
 *  ```
 *  Note that this function is called *after* the end_iteration_step
 *  handler.
 *
 * By default this class instructs the scf to run for a maximum of 100
 * iterations
 * and compute all eigenpairs.
 *
 * \see gscf::PlainScf
 */
template <typename ScfState>
struct PlainScfControl : public ScfControlBase<ScfState> {
  typedef ScfControlBase<ScfState> base_type;
  typedef typename base_type::size_type size_type;

  PlainScfControl()
        //          Number of eigenpairs,               maxiter
        : base_type{IterationConstants<size_type>::all, 100} {}
};

/** \name  Basic and plain SCF algorithm
 *
 * This just keeps running until maxiter is reached. If you want to add
 * a convergence criterion override the is_converged(state) function.
 *
 * \tparam ProblemMatrix  The type of the problem matrix to be solved
 * \tparam ScfState       The precise scf state type which is available
 *                        to is_converged and all handler functions.
 * \tparam ScfControl     The scf control type which contains parameters
 *                        as well as functionality to check convergence.
 */
template <typename ProblemMatrix,
          typename ScfState = PlainScfState<ProblemMatrix>,
          typename ScfControl = PlainScfControl<ScfState>>
class PlainScf : public ScfBase<ScfState, ScfControl> {
public:
  typedef ScfBase<ScfState, ScfControl> base_type;
  typedef typename base_type::scf_state_type scf_state_type;

  typedef typename scf_state_type::probmat_type probmat_type;
  typedef typename scf_state_type::scalar_type scalar_type;
  typedef typename scf_state_type::size_type size_type;
  typedef typename scf_state_type::matrix_type matrix_type;
  typedef typename scf_state_type::vector_type vector_type;

  static_assert(std::is_same<ProblemMatrix, probmat_type>::value,
                "The ProblemMatrix type specified and the one implicit in the "
                "ScfState class have to agree.");

  /** Run a Plain SCF eigenproblem with the provided problem matrix
   *  and overlap matrix and return the final state.
   *
   * \param assert_nofail If set to true, any exception indicating convergence
   * failure will be passed upwards, else it will be suppressed.
   * In this case the fail message of the returned state gives valuable
   * information why the iteration failed.
   *
   * \note The returned state contains information on whether the SCF
   * converged or failed.
   **/
  scf_state_type solve(probmat_type probmat_bb,
                       const matrix_type& overlapmat_bb,
                       bool assert_nofail = true) const override;
};

template <typename ProblemMatrix, typename ScfState, typename ScfControl>
typename PlainScf<ProblemMatrix, ScfState, ScfControl>::scf_state_type
PlainScf<ProblemMatrix, ScfState, ScfControl>::solve(
      probmat_type probmat_bb, const matrix_type& overlapmat_bb,
      bool assert_nofail) const {
  // Construct new state:
  scf_state_type state(std::move(probmat_bb), overlapmat_bb);

  try {
    // Iterate:
    while (!base_type::is_converged(state)) {
      base_type::start_iteration_step(state);
      base_type::solve_eigensystem(state);
      base_type::update_problem_matrix(state);
      base_type::end_iteration_step(state);
    }
  } catch (ExcScfFailedToConverge& e) {
    if (assert_nofail) throw;
    return state;
  }

  // Return the final state
  return state;
}

}  // gscf
