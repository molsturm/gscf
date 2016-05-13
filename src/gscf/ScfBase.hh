#pragma once
#include "IterationConstants.hh"
#include "ScfControlBase.hh"
#include "ScfStateBase.hh"
#include "eig_sym_hack.hh"
#include <linalgwrap/Exceptions.hh>
#include <type_traits>

namespace gscf {

DefException1(ExcScfFailedToConverge, char*, << "The SCF failed to converge: "
                                             << arg1);

/** Class to provide building block functions and a common interface
 * for SCF procedures
 *
 * \see gscf::ScfStateBase
 * \see gscf::ScfControlBase
 * */
template <typename ScfState>
class ScfBase {
  static_assert(IsScfState<ScfState>::value,
                "ScfState needs to be a type derived from ScfStateBase");

public:
  typedef ScfState scf_state_type;
  typedef typename scf_state_type::probmat_type probmat_type;
  typedef typename scf_state_type::scalar_type scalar_type;
  typedef typename scf_state_type::size_type size_type;
  typedef typename scf_state_type::matrix_type matrix_type;
  typedef typename scf_state_type::vector_type vector_type;

  /** \name Run an SCF
   */
  ///@{
  /**
   * \brief Run the actual SCF solver.
   *
   * Perform a sequence of SCF iterations and return as soon as the SCF
   * converged or an error occurred.
   *
   * The returned state gives access to the most recent state and
   * indicates whether the iteration was successful.
   *
   * \note Most people will find the is_converged() function
   * of the SCF state most helpful.
   */
  virtual scf_state_type solve(probmat_type probmat_bb,
                               const matrix_type& overlapmat_bb) const = 0;

  /** \brief Run the solver starting from an old SCF state.
   *
   * It is assumed, that the state is not failed.
   *
   * Implementing classes are advised to overload this method in order to
   * make use of the old state.
   */
  virtual scf_state_type solve(const scf_state_type& old_state) const;

  /** Call the solve() function of the scf solver. If the solve is
   * unsuccessful an instance of ExcScfNotConverged is thrown */
  scf_state_type solve_assert(probmat_type probmat_bb,
                              const matrix_type& overlapmat_bb) const;

  /** Call the solve() function of the scf solver. If the solve is
   * unsuccessful an instance of ExcScfNotConverged is thrown */
  scf_state_type solve_assert(const scf_state_type& old_state) const;
  ///@}

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
  ///@}

  /** Construct an ScfBase object from an ScfControlBase object, which
   * contains the iteration parameters */
  ScfBase(const ScfControlBase<scf_state_type>& scf_control);

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
  linalgwrap::SubscriptionPointer<const ScfControlBase<scf_state_type>>
        m_scf_control_ptr;
};

//
// ------------------------------------
//

template <typename ScfState>
inline typename ScfBase<ScfState>::scf_state_type
ScfBase<ScfState>::solve_assert(probmat_type probmat_bb,
                                const matrix_type& overlapmat_bb) const {
  scf_state_type state = solve(std::move(probmat_bb), overlapmat_bb);
  assert_throw(!state.is_failed(),
               ExcScfFailedToConverge(state.fail_reason().c_str()));
  return state;
}

template <typename ScfState>
inline typename ScfBase<ScfState>::scf_state_type ScfBase<ScfState>::solve(
      const scf_state_type& old_state) const {
  return solve(*old_state.problem_matrix_ptr(), old_state.overlap_matrix());
}

template <typename ScfState>
inline typename ScfBase<ScfState>::scf_state_type
ScfBase<ScfState>::solve_assert(const scf_state_type& old_state) const {
  scf_state_type state = solve(old_state);
  assert_throw(!state.is_failed(), ExcScfFailedToConverge(state.fail_reason()));
  return state;
}

template <typename ScfState>
inline ScfBase<ScfState>::ScfBase(
      const ScfControlBase<scf_state_type>& scf_control)
      : m_scf_control_ptr{
              linalgwrap::make_subscription(scf_control, "ScfBase")} {}

template <typename ScfState>
inline void ScfBase<ScfState>::start_iteration_step(scf_state_type& s) const {
  // Assert that SCF is not failed.
  assert_dbg(!s.is_failed(),
             linalgwrap::ExcInvalidState(
                   ("SCF has failed, reason: " + s.fail_reason()).c_str()));

  // Increase the iteration count:
  s.increase_iteration_count();

  // Call the pre step handler:
  before_iteration_step(s);
}

template <typename ScfState>
void ScfBase<ScfState>::solve_eigensystem(scf_state_type& s) const {

  // The number of eigenpairs to compute:
  size_type n_eigenpairs = m_scf_control_ptr->n_eigenpairs;
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
    // There was an error:
    s.fail("Inner eigenproblem failed to converge");
    return;
  }

  // Update state:
  s.eigenvectors_ptr() = new_eigenvectors_ptr;
  s.eigenvalues_ptr() = new_eigenvalues_ptr;

  // call the handler:
  on_update_eigenpairs(s);
}

template <typename ScfState>
void ScfBase<ScfState>::update_problem_matrix(scf_state_type& s) const {
  const std::string coeff_key = "evec_coefficients";

  // Store new coefficients in a parameter map:
  linalgwrap::ParameterMap m;
  m.update(coeff_key, s.eigenvectors_ptr());

  if (!s.problem_matrix_ptr().unique()) {
    // s is not the only thing referencing the object behind the
    // problem_matrix_ptr shared pointer, so we need to copy the object
    // first. This assuresthat objects that depend on the SCF state (i.e.
    // those which actually copyied the shared pointer to some internal
    // place) still access the old, non-updated problem matrix.
    //
    // A typical scenario for this is the Pulay-DIIS, where a linear
    // combination
    // of Fock matrices is used as a guess for the next problem matrix.
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

template <typename ScfState>
inline void ScfBase<ScfState>::end_iteration_step(scf_state_type& s) const {
  // Call the handler:
  after_iteration_step(s);

  if (s.n_iter_count() >= m_scf_control_ptr->max_iter) {
    s.fail("Maximum number of iterations reached.");
  }
}

}  // namespace gscf
