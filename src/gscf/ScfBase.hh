#pragma once
#include "ExcScfFailedToConverge.hh"
#include "IterationConstants.hh"
#include "ScfControlBase.hh"
#include "ScfStateBase.hh"
#include "eig_sym_hack.hh"
#include <type_traits>

namespace gscf {

/** Class to provide building block functions and a common interface
 * for SCF procedures
 *
 * \see gscf::ScfStateBase
 * \see gscf::ScfControlBase
 * */
template <typename ScfState, typename ScfControl>
class ScfBase {
  static_assert(IsScfState<ScfState>::value,
                "ScfState needs to be a type derived from ScfStateBase");

  static_assert(IsScfControl<ScfControl>::value,
                "ScfState needs to be a type derived from ScfControlBase");

  static_assert(
        std::is_base_of<typename ScfControl::scf_state_type, ScfState>::value,
        "The ScfState needs to be derived off ScfControl::scf_state_type for "
        "proper functioning");

public:
  typedef ScfControl scf_control_type;
  typedef ScfState scf_state_type;
  typedef typename scf_state_type::probmat_type probmat_type;
  typedef typename scf_state_type::scalar_type scalar_type;
  typedef typename scf_state_type::size_type size_type;
  typedef typename scf_state_type::matrix_type matrix_type;
  typedef typename scf_state_type::vector_type vector_type;

  /** \name Run an SCF
   */
  ///@{
  /** \brief Run the scf solver with the provided problem matrix
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
  virtual scf_state_type solve(probmat_type probmat_bb,
                               const matrix_type& overlapmat_bb,
                               bool assert_nofail) const = 0;

  /** \brief Run the solver starting from an old SCF state.
   *
   * It is assumed, that the input state is not failed.
   * Note that the fail bit can be unset using the clear_failed() function.
   *
   * Implementing classes are advised to overload this method in order to
   * make use of the old state.
   *
   * \param assert_nofail If set to true, any exception indicating convergence
   * failure will be passed upwards, else it will be suppressed.
   * In this case the fail message of the returned state gives valuable
   * information why the iteration failed.
   */
  virtual scf_state_type solve(const scf_state_type& old_state,
                               bool assert_nofail = true) const;
  ///@}

  /** Access scf control read-only */
  const scf_control_type scf_control() const;

  /** Access scf control read-write */
  scf_control_type scf_control();

protected:
  /** \name Handler functions
   * Various virtual handler functions, which are called when
   * certain events happen.
   */
  ///@{
  /* Handler which is called before an iteration step is performed
   *
   * The iteration count has already been incremented.
   * */
  virtual void before_iteration_step(scf_state_type&) const {}

  /** Handler which is called once an iteration step finishes
   *
   * This is the last thing called before the convergence and sanity
   * checks.
   * */
  virtual void after_iteration_step(scf_state_type&) const {}

  /** Handler which is called once the problem matrix has been diagonalised
   *  and new eigenpairs are obtained for this problem matrix
   */
  virtual void on_update_eigenpairs(scf_state_type&) const {}

  /** Handler which is called once the new problem matrix has been formed
   *  from the current set of eigenvectors
   */
  virtual void on_update_problem_matrix(scf_state_type&) const {}

  /** Handler which is called once the iteration has converged */
  virtual void on_converged(scf_state_type&) const {}

  /** Handler which is called by fail_scf once the iteration has failed */
  virtual void on_failed(scf_state_type&) const {}
  ///@}

  /** Construct an ScfBase object. */
  ScfBase() : m_scf_control{} {}

  /** \name Scf building blocks.
   * Various virtual handler functions, which are called when
   * certain events happen.
   */
  ///@{
  /** \brief Start the next iteration step
   *
   * Increase the iteration count and call the pre_step_handler.
   **/
  void start_iteration_step(scf_state_type& s) const;

  /** \brief Check whether convergence has been achieved
   * Return true if yes, else false
   **/
  bool is_converged(scf_state_type& s) const;

  /** \brief Fail an scf iteration.
   *
   * Set the fail message of the state to fail_message.
   * Then call the on_failed handler and throw an
   * ExcScfFailedToConverge with the same message.
   */
  void fail_scf(scf_state_type& s, const ScfFailReason& fail_reason) const;

  /** \brief Solve the eigensystem currently represented by the
   *  SCF state and place results back inside the state
   *
   *  After everything is done, we call the update_eigenpair_handler.
   *
   *  Notice, that we just re-assign the shared pointers
   *  and that their magic is in charge for cleaning up
   *  unused storage automatically.
   */
  void solve_eigensystem(scf_state_type& s) const;

  /** Update the problem matrix and replace the current
   *  one in the SCF state.
   *
   *  Notice, that we just re-assign the shared pointers
   *  and that their magic is in charge for cleaning up
   *  unused storage automatically.
   */
  void update_problem_matrix(scf_state_type& s) const;

  /** End the current iteration step.
   *
   *  In principle it check that we are not beyond
   *   max_iter
   */
  void end_iteration_step(scf_state_type& s) const;
  ///@}

private:
  //! The SCF control object, which controls the convergence
  //  checks and iterations done in this class.
  scf_control_type m_scf_control;
};

//
// ------------------------------------
//

template <typename ScfState, typename ScfControl>
inline typename ScfBase<ScfState, ScfControl>::scf_state_type
ScfBase<ScfState, ScfControl>::solve(const scf_state_type& old_state,
                                     bool assert_nofail) const {
  assert_dbg(!old_state.is_failed(),
             linalgwrap::ExcInvalidState("Cannot make use of a failed state"));
  return solve(*old_state.problem_matrix_ptr(), old_state.overlap_matrix(),
               assert_nofail);
}

template <typename ScfState, typename ScfControl>
const typename ScfBase<ScfState, ScfControl>::scf_control_type
ScfBase<ScfState, ScfControl>::scf_control() const {
  return m_scf_control;
}

template <typename ScfState, typename ScfControl>
typename ScfBase<ScfState, ScfControl>::scf_control_type
ScfBase<ScfState, ScfControl>::scf_control() {
  return m_scf_control;
}

template <typename ScfState, typename ScfControl>
inline void ScfBase<ScfState, ScfControl>::start_iteration_step(
      scf_state_type& s) const {
  // Assert that SCF is not failed.
  assert_dbg(!s.is_failed(),
             linalgwrap::ExcInvalidState(
                   ("SCF has failed, reason: " + s.fail_reason()).c_str()));

  // Increase the iteration count:
  s.increase_iteration_count();

  // Call the pre step handler:
  before_iteration_step(s);
}

template <typename ScfState, typename ScfControl>
bool ScfBase<ScfState, ScfControl>::is_converged(scf_state_type& s) const {
  bool converged = !s.is_failed() && m_scf_control.is_converged(s);
  if (converged) {
    on_converged(s);
  }
  return converged;
}

template <typename ScfState, typename ScfControl>
void ScfBase<ScfState, ScfControl>::fail_scf(
      scf_state_type& s, const ScfFailReason& fail_reason) const {
  s.fail(fail_reason.as_string());
  on_failed(s);
  assert_throw(false, ExcScfFailedToConverge(fail_reason));
}

template <typename ScfState, typename ScfControl>
void ScfBase<ScfState, ScfControl>::solve_eigensystem(scf_state_type& s) const {
  // Assert that SCF is not failed.
  assert_dbg(!s.is_failed(),
             linalgwrap::ExcInvalidState(
                   ("SCF has failed, reason: " + s.fail_reason()).c_str()));

  // The number of eigenpairs to compute:
  size_type n_eigenpairs = m_scf_control.n_eigenpairs;
  if (n_eigenpairs == IterationConstants<size_type>::all) {
    // Compute all eigenpairs:
    n_eigenpairs = s.problem_matrix_ptr()->n_cols();
  }

  assert_greater_equal(n_eigenpairs, s.problem_matrix_ptr()->n_cols());

  // Containers for the new eigenpairs:
  auto new_eigenvectors_ptr = std::make_shared<matrix_type>(
        s.problem_matrix_ptr()->n_rows(), n_eigenpairs);
  auto new_eigenvalues_ptr = std::make_shared<vector_type>(n_eigenpairs);

  // Do eigenproblem:
  // TODO Take the number of eigenpairs requested into account here
  bool res = detail::eig_sym_hack(*s.problem_matrix_ptr(), s.overlap_matrix(),
                                  *new_eigenvectors_ptr, *new_eigenvalues_ptr);

  if (!res) {
    fail_scf(s, ScfFailReason::ReasonId::INNER_EIGENSOLVER_FAILED);
    return;
  }

  // Update state:
  s.eigenvectors_ptr() = new_eigenvectors_ptr;
  s.eigenvalues_ptr() = new_eigenvalues_ptr;

  // call the handler:
  on_update_eigenpairs(s);
}

template <typename ScfState, typename ScfControl>
void ScfBase<ScfState, ScfControl>::update_problem_matrix(
      scf_state_type& s) const {
  // Assert that SCF is not failed.
  assert_dbg(!s.is_failed(),
             linalgwrap::ExcInvalidState(
                   ("SCF has failed, reason: " + s.fail_reason()).c_str()));

  // Store new coefficients in a parameter map:
  linalgwrap::ParameterMap m;
  m.update(m_scf_control.update_key, s.eigenvectors_ptr());

  if (!s.problem_matrix_ptr().unique()) {
    // s is not the only thing referencing the object behind the
    // problem_matrix_ptr shared pointer, so we need to copy the object
    // first. This assuresthat objects that depend on the SCF state (i.e.
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
          std::make_shared<probmat_type>(*s.problem_matrix_ptr());

    // Replace the current one in the state:
    s.problem_matrix_ptr() = new_problem_matrix_ptr;
  }

  // Update the new problem matrix:
  s.problem_matrix_ptr()->update(m);

  // call the handler:
  on_update_problem_matrix(s);
}

template <typename ScfState, typename ScfControl>
inline void ScfBase<ScfState, ScfControl>::end_iteration_step(
      scf_state_type& s) const {
  // Assert that SCF is not failed.
  assert_dbg(!s.is_failed(),
             linalgwrap::ExcInvalidState(
                   ("SCF has failed, reason: " + s.fail_reason()).c_str()));

  // Call the handler:
  after_iteration_step(s);

  if (s.n_iter_count() >= m_scf_control.max_iter) {
    fail_scf(s, ScfFailReason::ReasonId::MAXIMUM_ITERATIONS_REACHED);
  }
}

}  // namespace gscf
