#pragma once
#include "CircularBuffer.hh"
#include "ScfBase.hh"
#include "ScfControlBase.hh"
#include "ScfStateBase.hh"
#include "linsolve_hack.hh"

namespace gscf {

/** Default scf state for a Pulay DIIS SCF
 *
 * \tparam ProblemMatrix The type of the Problem matrix object.
 * */
template <typename ProblemMatrix>
// TODO invent type-traits to determine "linear combination type"
struct PulayDiisScfState
      : public ScfStateBase<ProblemMatrix,
                            linalgwrap::LazyMatrixSum<
                                  typename ProblemMatrix::stored_matrix_type>> {
  typedef ScfStateBase<
        ProblemMatrix,
        linalgwrap::LazyMatrixSum<typename ProblemMatrix::stored_matrix_type>>
        base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::diagmat_type diagmat_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  /** The last n_prev_steps eigenvector pointers,
   * front() is the oldest, back() the most recent*/
  CircularBuffer<std::shared_ptr<const matrix_type>> prev_eigenvectors_ptrs;

  /** The last n_prev_steps eigenvalue pointers,
   * front() is the oldest, back() the most recent*/
  CircularBuffer<std::shared_ptr<const vector_type>> prev_eigenvalues_ptrs;

  /** The last n_prev_steps problem matrix pointers,
   * front() is the oldest, back() the most recent */
  CircularBuffer<std::shared_ptr<const probmat_type>> prev_problem_matrix_ptrs;

  /** The last n_prev_steps error matrices
   * front() is the oldest, back() the most recent*/
  CircularBuffer<matrix_type> errors;

  /** Overlap between the errors of the previous steps
   * front() is the oldest set of overlaps, back() the most recent
   *
   * If front() points to the overlaps computed at the ith
   * step of the algorithm, then the buffer contains:
   * front    -> (   <i|i>     <i|i-1>   <i|i-2> ...   <i|i-n>   )
   * front+1  -> ( <i+1|i+1> <i+1|i>   <i+1|i-1> ... <i+1|i-n+1> )
   * ...
   * back     -> ( <i+n|i+n> <i+n|i+n-1> ...         <i+n|i>     )
   *
   * Of cause if i < n, then the first vector contains less then n
   * elements. Each iteration, only a new vector needs to be added at
   * the back of the buffer. No operation shall touch the elements
   * of any vector which are no longer valid (i.e. with <i|j> and
   * j < i).
   *
   * In all above <i|j> is the overlap between the error computed
   * at iteration step i and the error computed at iteration step j.
   **/
  CircularBuffer<std::vector<scalar_type>> error_overlaps;

  /** The current diis coefficients, in the same order as the problem matrices,
   * i.e.
   * the first coefficient belongs to the oldest problem matrix, the last
   * coefficient
   * to the most recent. */
  vector_type diis_coefficients;

  PulayDiisScfState(probmat_type probmat, const matrix_type& overlap_mat,
                    const size_type n_prev_steps)
        : base_type{std::move(probmat), overlap_mat},
          prev_eigenvectors_ptrs{n_prev_steps},
          prev_eigenvalues_ptrs{n_prev_steps},
          prev_problem_matrix_ptrs{n_prev_steps},
          errors{n_prev_steps},
          error_overlaps{n_prev_steps},
          diis_coefficients(0) {}
};

/** \brief Class to provide the basic set of parameters and the convergence
 *  checks the ScfBase class needs.
 *
 *  Consider overriding the is_converged function to provide a more sensible
 *  convergence check. This version just checks whether the maximal l2-norm
 *  of Frobenius norm (if the error is a vector or a matrix, respectively)
 *  does not reach beyond a certain threshold (by default 5*10^{-7}).
 *  The signature is
 *  ```
 *  bool is_converged(const scf_state_type&) const;
 *  ```
 *  Note that this function is called *after* the end_iteration_step
 *  handler.
 *
 * By default this class instructs the scf to run for a maximum of 100
 * iterations and compute all eigenpairs. For convergence the maximum a
 * diis error matrix element can have is 5*10^{-7}.
 *
 * \see gscf::PlainScf
 */
template <typename ScfState>
struct PulayDiisScfControl : public ScfControlBase<ScfState> {
  typedef ScfControlBase<ScfState> base_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::scf_state_type scf_state_type;
  typedef typename base_type::matrix_type matrix_type;

  //! The number of previous scf steps to consider
  size_type n_prev_steps;

  /** Maximum value the Frobenius norm of the
   *  most recent error vector/matrix may have. */
  scalar_type max_diis_error;

  bool is_converged(const scf_state_type& s) const override {
    if (s.errors.empty()) return false;

    const matrix_type& last_error = s.errors.back();
    return last_error.norm_frobenius() < max_diis_error;
  }

  PulayDiisScfControl()
        //          Number of eigenpairs,               maxiter
        : base_type{IterationConstants<size_type>::all, 100},
          n_prev_steps{5},
          max_diis_error(5e-7) {}
};

/** \name Pulay DIIS SCF algorithm
 *
 * Run the Pulay DIIS SCF algorithm until the convergence has been reached
 * (by default Frobenius norm of the error < 5e-7, see PulayDiisScfControl).
 *
 * By default the error is taken to be the residual error, i.e. the
 * difference between the subsequent eigenpairs calculated. Note, that this
 * is not a good error for the HF-SCF procedure! Override the calculate_error
 * function in order to provide a different error.
 *
 * By default the next DIIS guess is formed by linearly combining the
 *
 * \tparam ProblemMatrix  The type of the problem matrix to be solved
 * \tparam ScfState       The precise scf state type which is available
 *                        to is_converged and all handler functions.
 * \tparam ScfControl     The scf control type which contains parameters
 *                        as well as functionality to check convergence.
 */
template <typename ProblemMatrix,
          typename ScfState = PulayDiisScfState<ProblemMatrix>,
          typename ScfControl = PulayDiisScfControl<ScfState>>
class PulayDiisScf : public ScfBase<ScfState, ScfControl> {
public:
  typedef ScfBase<ScfState, ScfControl> base_type;
  typedef typename base_type::scf_state_type scf_state_type;

  typedef typename scf_state_type::probmat_type probmat_type;
  typedef typename scf_state_type::diagmat_type diagmat_type;
  typedef typename scf_state_type::scalar_type scalar_type;
  typedef typename scf_state_type::size_type size_type;
  typedef typename scf_state_type::matrix_type matrix_type;
  typedef typename scf_state_type::vector_type vector_type;

  static_assert(std::is_same<ProblemMatrix, probmat_type>::value,
                "The ProblemMatrix type specified and the one implicit in the "
                "ScfState class have to agree.");

  static_assert(
        std::is_base_of<PulayDiisScfState<ProblemMatrix>, ScfState>::value,
        "ScfState needs to be derived off PulayDiisScfState");
  static_assert(
        std::is_base_of<PulayDiisScfControl<ScfState>, ScfControl>::value,
        "ScfControl needs to be derived off PulayDiisScfControl");

public:
  /** Run a Pulay DIIS SCF eigenproblem with the provided problem matrix
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

protected:
  /** \name Handler functions
   * Various virtual handler functions, which are called when
   * certain events happen.
   */
  ///@{
  /* Handler which is called once the new guess for the self-consistent
   * problem matrix has been formed by the DIIS.
   *
   * \note The guess can be obtained via state.diagonalised_matrix_ptr().
   * */
  virtual void on_new_diis_diagmat(scf_state_type&) const {}
  ///@}

  /** Calculate the error from a provided scf state, which should be
   *  a PulayDiisScfState.
   *
   *  This implementation calculates the residual error, i.e. the difference
   *  vectors between the current and the most recent eigenvectors.
   *
   * \note This function is called once the updated problem matrix has
   * been obtained.
   **/
  virtual matrix_type calculate_error(const scf_state_type& s) const;

  /** Calculate the overlap between the most recent error matrix and
   * all previous error matrices, which are stored.
   * The values returned are in the order <i|i>, <i|i-1>, ... <i|i-n>
   *
   * \note This function is called after the new error for this iteration
   * has been computed.
   */
  void append_new_overlaps(scf_state_type& s) const;

  /** Use the error overlaps of the current state to calculate the
   *  new diis coefficients and store them inside the state.
   *
   * \note This function is called as the first thing
   * once the iteration has started
   */
  void update_diis_coefficients(scf_state_type& s) const;

  /** Use the current diis coefficient in order to compute the new
   *  DIIS guess for the self-consistent problem matrix.
   *  The result is stored inside diagonalised_matrix_ptr
   *
   *  Then on_new_diis_diagmat is called.
   *
   *  \note This function is called once the new diis coefficients
   *  have been obtained and before diagonalisation.
   */
  void update_diis_diagmat(scf_state_type& s) const;

  /** Compute the system matrix for the DIIS linear system to be solved.
   * This matrix consists of the error overlaps in the upper-left corner
   * and a row of -1 and a column of -1 following thereafter. It looks like
   * \f[
   *    \left(
   *    begin{array}{ccccc}
   *    B_{00} & B_{01} & \cdots & B_{0,m-1} & -1 \\
   *    B_{10} & B_{11} & \cdots & B_{1,m-1} & -1 \\
   *    \vdots & \vdots & \vdots & \vdots    & \vdots \\
   *    B_{m-1,0} & B_{m-1,1} & \cdots & B_{m-1,m-1} & -1 \\
   *    -1     &  -1    & \cdots &  -1       &  0
   *    end{array}
   *    right)
   * \f]
   * where \f$B_{00}\f$ is the overlap between error matrices 0 and 0.
   * The order is from oldest to most recent, i.e. \f$B_{00}\f$ refers
   * to the oldest problem matrix and \f$B_{m-1,m-1}\f$ to the most recent.
   * */
  matrix_type diis_linear_system_matrix(const scf_state_type& s) const;
};

//
// ------------------------------------------
//

template <typename ProblemMatrix, typename ScfState, typename ScfControl>
typename PulayDiisScf<ProblemMatrix, ScfState, ScfControl>::scf_state_type
PulayDiisScf<ProblemMatrix, ScfState, ScfControl>::solve(
      probmat_type probmat_bb, const matrix_type& overlapmat_bb,
      bool assert_nofail) const {
  // Construct new state:
  scf_state_type state(std::move(probmat_bb), overlapmat_bb,
                       base_type::scf_control().n_prev_steps);

  try {
    while (!base_type::is_converged(state)) {
      base_type::start_iteration_step(state);

      // Compute new DIIS coefficients and guess
      update_diis_coefficients(state);
      update_diis_diagmat(state);

      // Solve the eigensystem of the DIIS guess and push back results
      base_type::update_eigenpairs(state);
      state.prev_eigenvectors_ptrs.push_back(state.eigenvectors_ptr());
      state.prev_eigenvalues_ptrs.push_back(state.eigenvalues_ptr());

      // The DIIS guess may contain copies of all operators in history.
      // To release memory, we reset it at this stage
      state.diagonalised_matrix_ptr().reset();

      // Form the new problem matrix and append to the history
      base_type::update_problem_matrix(state);
      state.prev_problem_matrix_ptrs.push_back(state.problem_matrix_ptr());

      // Compute new errors and error overlaps.
      state.errors.push_back(calculate_error(state));
      append_new_overlaps(state);

      base_type::end_iteration_step(state);
    }
  } catch (ExcScfFailedToConverge& e) {
    if (assert_nofail) throw;
    return state;
  }

  // Return the final state
  return state;
}

template <typename ProblemMatrix, typename ScfState, typename ScfControl>
typename PulayDiisScf<ProblemMatrix, ScfState, ScfControl>::matrix_type
PulayDiisScf<ProblemMatrix, ScfState, ScfControl>::diis_linear_system_matrix(
      const scf_state_type& s) const {
  size_type n_errors = s.error_overlaps.size();

  // The B Matrix of the DIIS (system matrix of the linear system to be solved)
  matrix_type B(n_errors + 1, n_errors + 1, false);

  auto it = std::begin(s.error_overlaps);
  for (size_type i = 0; i < n_errors; ++i, ++it) {
    // The diagonal element <i|i>
    B(i, i) = (*it)[0];

    for (size_type j = 1; j < i + 1; ++j) {
      // The off-diagonals => fill in symmetric manor
      B(i, i - j) = B(i - j, i) = (*it)[j];
    }

    // The extra elements of the B matrix (to fulfil the constraint
    // that the sum of all DIIS coefficients should be 1)
    B(n_errors, i) = B(i, n_errors) = -1.;
  }

  // Fill the remaining B matrix element in the extra row.
  B(n_errors, n_errors) = 0.;

  return B;
}

template <typename ProblemMatrix, typename ScfState, typename ScfControl>
void PulayDiisScf<ProblemMatrix, ScfState,
                  ScfControl>::update_diis_coefficients(scf_state_type& s)
      const {
  assert_dbg(!s.is_failed(),
             linalgwrap::ExcInvalidState(
                   ("SCF has failed, reason: " + s.fail_reason()).c_str()));

  size_type n_errors = s.error_overlaps.size();

  if (n_errors == 0) {
    // Store empty vector of coefficients and return.
    s.diis_coefficients = vector_type(0);
    return;
  }
  if (n_errors == 1) {
    // Store a single 1 as coefficient and return.
    s.diis_coefficients = vector_type{1.};
    return;
  }

  // The B Matrix of the DIIS (system matrix of the linear system to be solved)
  matrix_type B = diis_linear_system_matrix(s);

  // Build the RHS
  vector_type rhs(n_errors + 1);  // initialise with zeros
  rhs(n_errors) = -1;

  // Do the linear solve
  vector_type x(n_errors + 1);
  bool success = detail::linsolve_hack(B, rhs, x);
  if (!success) {
    base_type::fail_scf(s, ScfFailReason::ReasonId::DIIS_STEP_FAILED);
  }

  // Keep the first n_errors entries of the vector.
  // TODO more clever way to do this?
  s.diis_coefficients =
        vector_type(std::begin(x), std::next(std::begin(x), n_errors));
}

template <typename ProblemMatrix, typename ScfState, typename ScfControl>
void PulayDiisScf<ProblemMatrix, ScfState, ScfControl>::update_diis_diagmat(
      scf_state_type& s) const {
  using namespace linalgwrap;

  assert_dbg(!s.is_failed(),
             linalgwrap::ExcInvalidState(
                   ("SCF has failed, reason: " + s.fail_reason()).c_str()));
  assert_dbg(s.prev_problem_matrix_ptrs.size() == s.diis_coefficients.size(),
             linalgwrap::ExcInternalError());

  // Initialise an empty matrix for the DIIS guess in the state:
  s.diagonalised_matrix_ptr() = std::make_shared<diagmat_type>();

  // Build DIIS guess matrix to be diagonalised:
  diagmat_type& new_diis_guess = *s.diagonalised_matrix_ptr();
  if (s.diis_coefficients.size() == 0) {
    assert_dbg(s.error_overlaps.size() == 0, ExcInternalError());
    // This is the first run and there are no coefficients:
    new_diis_guess += *s.problem_matrix_ptr();
  } else {
    // Form linear combination according to coefficients:
    auto probmat_pit = std::begin(s.prev_problem_matrix_ptrs);
    size_t i = 0;
    for (; probmat_pit != std::end(s.prev_problem_matrix_ptrs);
         ++probmat_pit, ++i) {
      const probmat_type& mat = **probmat_pit;
      new_diis_guess += s.diis_coefficients[i] * mat;
    }
    assert_dbg(i == s.diis_coefficients.size(), ExcInternalError());
  }

  on_new_diis_diagmat(s);
}

template <typename ProblemMatrix, typename ScfState, typename ScfControl>
typename PulayDiisScf<ProblemMatrix, ScfState, ScfControl>::matrix_type
PulayDiisScf<ProblemMatrix, ScfState, ScfControl>::calculate_error(
      const scf_state_type& s) const {
  assert_dbg(!s.is_failed(),
             linalgwrap::ExcInvalidState(
                   ("SCF has failed, reason: " + s.fail_reason()).c_str()));

  const matrix_type& prev_evec = *s.prev_eigenvectors_ptrs.back();
  const matrix_type& cur_evec = *s.eigenvectors_ptr();
  return cur_evec - prev_evec;
}

template <typename ProblemMatrix, typename ScfState, typename ScfControl>
void PulayDiisScf<ProblemMatrix, ScfState, ScfControl>::append_new_overlaps(
      scf_state_type& s) const {
  assert_dbg(!s.is_failed(),
             linalgwrap::ExcInvalidState(
                   ("SCF has failed, reason: " + s.fail_reason()).c_str()));

  // The errors in s.errors are ordered from oldest(front) to newest(back)
  // In the overlap vector we need the order
  //      <n|n>, <n|n-1>, <n|n-m+1>, ..., <n|n-m>
  // where n is the current iteration number.

  // Resulting overlaps:
  std::vector<scalar_type> ret(s.errors.size(),
                               linalgwrap::Constants<scalar_type>::zero);
  // Fill vector from behind (reverse iterator)
  auto itoverlap = std::rbegin(ret);

  const matrix_type& cur_error = s.errors.back();
  for (const auto& error : s.errors) {
    // go through errors from front(old) to back(new)
    // and fill ret vector from behind.

    assert_size(error.n_rows(), cur_error.n_rows());
    assert_size(error.n_cols(), cur_error.n_cols());

    // TODO There really should be a function in linalgwrap for that.
    // like a dot_with between matrices or a function to view a matrix
    // as a linear vector and then have a dot_with on that
    for (auto it = std::begin(cur_error); it != std::end(cur_error); ++it) {
      *itoverlap += (*it) * error(it.row(), it.col());
    }
    assert_finite(*itoverlap);
    ++itoverlap;
  }

  s.error_overlaps.push_back(ret);

  std::cout << "   overlap:" << std::endl;
  for (auto it = std::begin(s.error_overlaps); it != std::end(s.error_overlaps);
       ++it) {
    std::cout << "    ";
    std::copy(std::begin(*it), std::end(*it),
              std::ostream_iterator<scalar_type>(std::cout, " "));
    std::cout << std::endl;
  }
}

}  // namespace gscf
