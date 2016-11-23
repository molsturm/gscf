#pragma once
#include "PulayDiisScfKeys.hh"
#include "ScfBase.hh"
#include "ScfStateBase.hh"
#include <krims/CircularBuffer.hh>
#include <linalgwrap/solve.hh>

namespace gscf {

DefSolverException1(ExcDiisStepFailed, std::string, details,
                    << "The SCF DIIS convergence accelerator step has failed: "
                    << details);

/** Default scf state for a Pulay DIIS SCF
 *
 * \tparam ProblemMatrix The type of the Problem matrix object.
 * */
template <typename ProblemMatrix, typename OverlapMatrix>
// TODO invent type-traits to determine "linear combination type"
struct PulayDiisScfState
      : public ScfStateBase<ProblemMatrix, OverlapMatrix,
                            linalgwrap::LazyMatrixSum<
                                  typename ProblemMatrix::stored_matrix_type>> {
  typedef ScfStateBase<
        ProblemMatrix, OverlapMatrix,
        linalgwrap::LazyMatrixSum<typename ProblemMatrix::stored_matrix_type>>
        base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::diagmat_type diagmat_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  /** The last n_prev_steps eigenvector pointers,
   * front() is the oldest, back() the most recent*/
  krims::CircularBuffer<
        std::shared_ptr<const linalgwrap::MultiVector<vector_type>>>
        prev_eigenvectors_ptrs;

  /** The last n_prev_steps eigenvalue pointers,
   * front() is the oldest, back() the most recent*/
  krims::CircularBuffer<std::shared_ptr<const std::vector<scalar_type>>>
        prev_eigenvalues_ptrs;

  /** The last n_prev_steps problem matrix pointers,
   * front() is the oldest, back() the most recent */
  krims::CircularBuffer<std::shared_ptr<const probmat_type>>
        prev_problem_matrix_ptrs;

  /** The last n_prev_steps error matrices
   * front() is the oldest, back() the most recent*/
  krims::CircularBuffer<matrix_type> errors;

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
  krims::CircularBuffer<std::vector<scalar_type>> error_overlaps;

  /** The current diis coefficients, in the same order as the problem matrices,
   * i.e.
   * the first coefficient belongs to the oldest problem matrix, the last
   * coefficient
   * to the most recent. */
  vector_type diis_coefficients;

  /** Construct an empty state, where all Circular buffers hold no elements */
  PulayDiisScfState(probmat_type probmat, const overlap_type& overlap_mat)
        : base_type{std::move(probmat), overlap_mat},
          prev_eigenvectors_ptrs{0},
          prev_eigenvalues_ptrs{0},
          prev_problem_matrix_ptrs{0},
          errors{0},
          error_overlaps{0},
          diis_coefficients(0) {}

  /** Resize all circular buffers to hold up to the given number
   * of previous SCF steps */
  void resize_buffers(const size_type n_prev_steps) {
    prev_eigenvectors_ptrs.max_size(n_prev_steps);
    prev_eigenvalues_ptrs.max_size(n_prev_steps);
    prev_problem_matrix_ptrs.max_size(n_prev_steps);
    errors.max_size(n_prev_steps);
    error_overlaps.max_size(n_prev_steps);
  }
};

/** \name Pulay DIIS SCF algorithm
 *
 * Run the Pulay DIIS SCF algorithm until the convergence has been reached
 * (by default Frobenius norm of the error < 5e-7, see PulayDiisScfControl).
 *
 * By default the error is taken to be the residual error, i.e. the
 * difference between the subsequent eigenpairs calculated. Note, that this
 * is not a good error for the HF-SCF procedure! Override the calculate_error
 * function in order to provide a different error estimation function.
 *
 * By default this error threshold is also used for convergence check
 * (this can be changed by overriding
 *  ```
 *  bool is_converged(const scf_state_type&) const;
 *  ```
 *  Note that this function is called *after* the end_iteration_step
 *  handler.
 *
 * By default the next DIIS guess is formed by linearly combining the last
 * 5 Fock matrices. See below how to change this.
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
 *   - n_prev_steps: The number of steps to consider in the DIIS.
 *                   (Default: 5)
 *   - max_error_norm:  If the Frobenius norm of the most recent error
 *                      vector/matrix computed by ``calculate_error``
 *                      is below this value, we consider the iteration
 *                      converged.
 *
 * \tparam ScfState  The precise scf state type which is available
 *                   to is_converged and all handler functions.
 *                   This also fixes the type of the problem matrix
 *                   and of the overlap matrix. Should be a the class
 *                   PulayDiisScfState or a subclass of it.
 */
template <typename ScfState>
class PulayDiisScf : public ScfBase<ScfState> {
public:
  typedef ScfBase<ScfState> base_type;
  typedef typename base_type::state_type state_type;

  typedef typename state_type::probmat_type probmat_type;
  typedef typename state_type::overlap_type overlap_type;
  typedef typename state_type::diagmat_type diagmat_type;
  typedef typename state_type::scalar_type scalar_type;
  typedef typename state_type::real_type real_type;
  typedef typename state_type::size_type size_type;
  typedef typename state_type::matrix_type matrix_type;
  typedef typename state_type::vector_type vector_type;

  static_assert(std::is_base_of<PulayDiisScfState<probmat_type, overlap_type>,
                                ScfState>::value,
                "ScfState needs to be derived off PulayDiisScfState");

public:
  /** \name Constructor */
  //@{
  /** Construct a DIIS SCF solver with default parameters */
  PulayDiisScf() {}

  /** Construct a DIIS SCF solver setting the parameters from the map */
  PulayDiisScf(const krims::ParameterMap& map) : PulayDiisScf() {
    update_control_params(map);
  }
  //@}

  /** \name Iteration control */
  ///@{
  //! The number of previous scf steps to consider
  size_type n_prev_steps = 5;

  /** Maximum value the Frobenius norm of the
   *  most recent error vector/matrix may have. */
  real_type max_error_norm = 5e-7;

  bool is_converged(const state_type& s) const override {
    if (s.errors.empty()) return false;

    // Norm of the last error (s.errors.back) is below convergence
    // threshold
    return norm_frobenius(s.errors.back()) < max_error_norm;
  }

  /** Update control parameters from Parameter map */
  void update_control_params(const krims::ParameterMap& map) {
    base_type::update_control_params(map);
    n_prev_steps = map.at(PulayDiisScfKeys::n_prev_steps, n_prev_steps);
    max_error_norm = map.at(PulayDiisScfKeys::max_error_norm, max_error_norm);
  }
  ///@}

  /** Implementation of the SolverBase method */
  void solve_state(state_type& state) const override;

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
  virtual void on_new_diis_diagmat(state_type&) const {}
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
  virtual matrix_type calculate_error(const state_type& s) const;

  /** Calculate the overlap between the most recent error matrix and
   * all previous error matrices, which are stored.
   * The values returned are in the order <i|i>, <i|i-1>, ... <i|i-n>
   *
   * \note This function is called after the new error for this iteration
   * has been computed.
   */
  void append_new_overlaps(state_type& s) const;

  /** Use the error overlaps of the current state to calculate the
   *  new diis coefficients and store them inside the state.
   *
   * \note This function is called as the first thing
   * once the iteration has started
   */
  void update_diis_coefficients(state_type& s) const;

  /** Use the current diis coefficient in order to compute the new
   *  DIIS guess for the self-consistent problem matrix.
   *  The result is stored inside diagonalised_matrix_ptr
   *
   *  Then on_new_diis_diagmat is called.
   *
   *  \note This function is called once the new diis coefficients
   *  have been obtained and before diagonalisation.
   */
  void update_diis_diagmat(state_type& s) const;

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
  matrix_type diis_linear_system_matrix(const state_type& s) const;
};

//
// ------------------------------------------
//

template <typename ScfState>
void PulayDiisScf<ScfState>::solve_state(state_type& state) const {
  // Resize the buffers in the state to hold
  // the appropriate number of recent SCF steps
  state.resize_buffers(n_prev_steps);

  while (!base_type::convergence_reached(state)) {
    base_type::start_iteration_step(state);

    // Compute new DIIS coefficients and guess
    update_diis_coefficients(state);
    update_diis_diagmat(state);

    // Solve the eigensystem of the DIIS guess and push back results
    base_type::update_eigenpairs(state);
    state.prev_eigenvectors_ptrs.push_back(state.eigenvectors_ptr());
    state.prev_eigenvalues_ptrs.push_back(state.eigenvalues_ptr());

    // The DIIS guess (stored in the diagonalised_matrix_ptr)
    // may contain copies of all operators in history.
    // To release the memory of what we don't need any more,
    // we reset it at this stage
    state.diagonalised_matrix_ptr().reset();

    // Form the new problem matrix and append to the history
    base_type::update_problem_matrix(state);
    state.prev_problem_matrix_ptrs.push_back(state.problem_matrix_ptr());

    // Compute new errors and error overlaps.
    state.errors.push_back(calculate_error(state));
    append_new_overlaps(state);

    base_type::end_iteration_step(state);
  }
}

template <typename ScfState>
typename PulayDiisScf<ScfState>::matrix_type
PulayDiisScf<ScfState>::diis_linear_system_matrix(const state_type& s) const {
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

template <typename ScfState>
void PulayDiisScf<ScfState>::update_diis_coefficients(state_type& s) const {
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
  vector_type rhs(n_errors + 1, true);  // initialise with zeros
  rhs(n_errors) = -1;

  try {
    // Do the linear solve

    // TODO: Here one could·specify a varying tolerance depending on
    // how accurate we want to have the final result and depending on
    // how large our SCF error currently is. This is just a shot
    // and currently does nothing.
    krims::ParameterMap params{{"tolerance", max_error_norm / 100.}};

    vector_type x(n_errors + 1, false);  // no initialisation
    linalgwrap::solve_hermitian(B, x, rhs, params);

    // Keep the first n_errors entries of the vector.
    // TODO more clever way to do this?
    s.diis_coefficients =
          vector_type(std::begin(x), std::next(std::begin(x), n_errors));
  } catch (const linalgwrap::SolverException& e) {
    std::stringstream ss;
    e.print_extra(ss);
    solver_assert(false, s, ExcDiisStepFailed(ss.str()));
  }
}

template <typename ScfState>
void PulayDiisScf<ScfState>::update_diis_diagmat(state_type& s) const {
  using namespace linalgwrap;
  using namespace krims;
  assert_dbg(s.prev_problem_matrix_ptrs.size() == s.diis_coefficients.size(),
             ExcInternalError());

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

template <typename ScfState>
typename PulayDiisScf<ScfState>::matrix_type
PulayDiisScf<ScfState>::calculate_error(const state_type& s) const {
  typedef linalgwrap::MultiVector<vector_type> mvec_type;
  const mvec_type& prev_evec = *s.prev_eigenvectors_ptrs.back();
  const mvec_type& cur_evec = *s.eigenvectors_ptr();

  // TODO until something better exists in linalgwrap
  // (like setting individual columns or so)
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
void PulayDiisScf<ScfState>::append_new_overlaps(state_type& s) const {
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

    *itoverlap = dot(cur_error, error);
    assert_finite(*itoverlap);
    ++itoverlap;
  }

  s.error_overlaps.push_back(ret);
}

}  // namespace gscf
