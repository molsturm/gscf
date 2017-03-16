#pragma once
#include "FocklikeMatrix_i.hh"
#include "ScfBase.hh"
#include "ScfStateBase.hh"
#include "TruncatedOptDampScfKeys.hh"
#include <krims/CircularBuffer.hh>

namespace gscf {

/** Default scf state for a truncated optimal damping algorithm SCF
 *
 * \note Do not rely on the objects in this class. This is currently only
 *       a test implementation which can only do a fixed n_prev_steps ==2.
 *       The interface is very likely going to change once a generalised
 *       ODA is going to be implemented.
 *
 * \tparam ProblemMatrix The type of the Problem matrix object.
 * */
template <typename ProblemMatrix, typename OverlapMatrix>
struct TruncatedOptDampScfState
      : public ScfStateBase<
              ProblemMatrix, OverlapMatrix,
              linalgwrap::LazyMatrixSum<typename ProblemMatrix::stored_matrix_type>> {

  static_assert(IsFocklikeMatrix<ProblemMatrix>::value,
                "The ProblemMatrix has to be a fock-like matrix for this SCF algorithm.");

  typedef ScfStateBase<
        ProblemMatrix, OverlapMatrix,
        linalgwrap::LazyMatrixSum<typename ProblemMatrix::stored_matrix_type>>
        base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::diagmat_type diagmat_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;
  typedef typename base_type::esoln_type esoln_type;

  /** The problem matrix of the previous step */
  std::shared_ptr<const probmat_type> prev_problem_matrix_ptr;

  /** The determined damping coefficient to combine the
   *  previous and the current problem matrix giving the
   *  matrix to be diagonalised.
   *
   *  \f[ F^\text{diag} = (1-c_\text{damp}) F^\text{prev} + c_\text{damp} F^\text{curr}
   * \f]
   */
  scalar_type damping_coeff;

  /** Construct an empty state */
  TruncatedOptDampScfState(probmat_type probmat, const overlap_type& overlap_mat)
        : base_type{std::move(probmat), overlap_mat},
          prev_problem_matrix_ptr{nullptr},
          damping_coeff{1} {}
};

/** \name Truncated optimal damping SCF algorithm
 *
 * TODO update this section
 *
 * Run the Pulay DIIS SCF algorithm until the convergence has been reached
 * (by default Frobenius norm of the error < 5e-7).
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
 *                      converged (default: 5e-7)
 *
 * \tparam ScfState  The precise scf state type which is available
 *                   to is_converged and all handler functions.
 *                   This also fixes the type of the problem matrix
 *                   and of the overlap matrix. Should be a the class
 *                   TruncatedOptDampScfState or a subclass of it.
 */
template <typename ScfState>
class TruncatedOptDampScf : public ScfBase<ScfState> {
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

  static_assert(std::is_base_of<TruncatedOptDampScfState<probmat_type, overlap_type>,
                                ScfState>::value,
                "ScfState needs to be derived off TruncatedOptDampScfState");

 public:
  /** \name Constructor */
  //@{
  /** Construct a DIIS SCF solver with default parameters */
  TruncatedOptDampScf() {}

  /** Construct a DIIS SCF solver setting the parameters from the map */
  TruncatedOptDampScf(const krims::GenMap& map) : TruncatedOptDampScf() {
    update_control_params(map);
  }
  //@}

  /** \name Iteration control */
  ///@{
  /** The number of previous scf steps to consider in the
   *  expansion to determine the optimal damping factor */
  size_type n_prev_steps = 2;

  /** Update control parameters from Parameter map */
  void update_control_params(const krims::GenMap& map) {
    base_type::update_control_params(map);
    n_prev_steps = map.at(TruncatedOptDampScfKeys::n_prev_steps, n_prev_steps);
  }

  /** Get the current settings of all internal control parameters and
   *  update the GenMap accordingly.
   */
  void get_control_params(krims::GenMap& map) const {
    base_type::get_control_params(map);
    map.update(TruncatedOptDampScfKeys::n_prev_steps, n_prev_steps);
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
   * problem matrix has been formed by the ODA algorithm.
   *
   * \note The guess can be obtained via state.diagonalised_matrix().
   * */
  virtual void on_new_diagmat(state_type&) const {}
  ///@}

  /** Use the problem-matrix-coefficient trace values of the current
   *  state to calculate the new damping coefficients.
   *
   * \note This function is called as the first thing
   * once the iteration has started
   */
  void update_damping_coefficients(state_type& s) const;

  /** Use the current damping coefficients in order to compute the new
   *  guess for the self-consistent problem matrix.
   *  The result is stored inside diagonalised_matrix_ptr
   *
   *  Then on_new_diagmat is called.
   *
   *  \note This function is called once the new damping coefficients
   *  have been obtained and before diagonalisation.
   */
  void update_oda_diagmat(state_type& s) const;
};

//
// ---------------------------------------------------------------
//

template <typename ScfState>
void TruncatedOptDampScf<ScfState>::solve_state(state_type& state) const {
  assert_dbg(!state.is_failed(), krims::ExcInvalidState("Cannot solve a failed state"));

  assert_dbg(n_prev_steps == 2, krims::ExcNotImplemented());

  while (!base_type::convergence_reached(state)) {
    base_type::start_iteration_step(state);

    // Compute new damping coefficients
    update_damping_coefficients(state);
    update_oda_diagmat(state);

    // Solve the eigensystem of the tODA problem matrix
    base_type::update_eigenpairs(state);

    // Update the matrix applies, noting that there are
    // various terms in the diagmat matrix we actually diagonalise
    const size_t n_terms = state.prev_problem_matrix_ptr == nullptr ? 1ul : 2ul;
    state.n_mtx_applies() += n_terms * state.eigenproblem_stats().n_mtx_applies();

    // The ODA guess (stored in the diagonalised_matrix_ptr)
    // may contain copies of all operators in history.
    // To release the memory of what we don't need any more,
    // we reset it at this stage
    state.diagonalised_matrix_ptr.reset();

    // Purge the problem matrix with the smallest
    // coefficient (lambda or (1-lambda) and update
    // the damping coefficient accordingly.
    if (state.damping_coeff > 0.5) {
      state.prev_problem_matrix_ptr = state.problem_matrix_ptr;
    } else {
      state.damping_coeff = 1 - state.damping_coeff;
    }

    // Form the new problem matrix
    base_type::update_problem_matrix(state);

    // Calculate the new error.
    base_type::update_last_error_norm(state);

    base_type::end_iteration_step(state);
  }
}

template <typename ScfState>
void TruncatedOptDampScf<ScfState>::update_damping_coefficients(state_type& s) const {
  if (s.prev_problem_matrix_ptr == nullptr) {
    s.damping_coeff = 1;
    return;
  }

  // The value of the trace $\tr F^{(n)} P^{(n+1)}$
  const scalar_type tr_cpc = [&s]() {
    // Aliase
    const auto& coeff_bf = s.eigensolution().evectors();
    const auto& fock_bb = s.problem_matrix();
    const auto& orben_f = s.eigensolution().evalues();

    // The range of occupied alphas
    const krims::Range<size_t> occ_a =
          s.problem_matrix().indices_subspace(OrbitalSpace::OCC_ALPHA);
    const krims::Range<size_t> occ_b =
          s.problem_matrix().indices_subspace(OrbitalSpace::OCC_BETA);

    // checks
    assert_greater(occ_a.length(), orben_f.size());
    assert_greater(occ_b.length(), orben_f.size());
    assert_sufficiently_tested(occ_a.length() == occ_b.length());

    // In the special case where the damping coefficient is 1 the diagonalised
    // matrix equals the previous problem matrix exactly. Therefore we do not need to
    // perform the trace in this special case and can just sum the eigenvalues of the
    // occupied orbitals to get the trace $\tr F^{(n)} P^{(n+1)}$
    if (s.damping_coeff == 1) {
      // Compute sum of the occupied alpha orbital energies
      const scalar_type sum_a =
            std::accumulate(std::begin(orben_f) + occ_a.lower_bound(),
                            std::begin(orben_f) + occ_a.upper_bound(), scalar_type(0));

      // Either return twice this value or add the occupied beta orbital energies
      return occ_a == occ_b
                   ? 2. * sum_a
                   : std::accumulate(std::begin(orben_f) + occ_b.lower_bound(),
                                     std::begin(orben_f) + occ_b.upper_bound(), sum_a);
    } else {
      // TODO This is poor-mans multiplexing. I feel this should go once we have
      // block-diagonality in place.
      // TODO increase apply count

      // Occupied coefficients (alpha)
      const auto ca_bo = coeff_bf.subview(occ_a);

      // Apply fock to occupied (alpha) orbitals
      //     -- O(n_bas*n_bas*n_occ)
      const auto Fca_bo = fock_bb * ca_bo;

      // Form the outer product sum and trace it
      //  -- O(n_bas*n_bas*n_occ)
      const auto tr_aa = trace(outer_prod_sum(ca_bo, Fca_bo));

      if (occ_a == occ_b) {
        return 2. * tr_aa;
      } else {
        // Repeat the same in the beta blocks

        const auto cb_bo = coeff_bf.subview(occ_b);
        const auto Fcb_bo = fock_bb * coeff_bf.subview(occ_b);
        const auto tr_bb = trace(outer_prod_sum(cb_bo, Fcb_bo));
        return tr_aa + tr_bb;
      }
    }
  }();

  // Compute s and c for the 2-step truncated ODA SCF
  const scalar_type oda_s = tr_cpc - s.prev_problem_matrix_ptr->energy_1e_terms() -
                            2. * s.prev_problem_matrix_ptr->energy_2e_terms();
  const scalar_type oda_c = s.problem_matrix().energy_1e_terms() - tr_cpc +
                            s.problem_matrix().energy_2e_terms() +
                            s.prev_problem_matrix_ptr->energy_2e_terms();

  // If the error in the SCF gets too small the numerics makes the values of oda_s and
  // oda_c go mad.
  // This catches the problematic cases and disables the ODA.
  if (oda_s >= 0 || fabs(oda_s) < 1e-14 || fabs(oda_c) < 1e-14) {
    s.damping_coeff = 1;
    return;
  }

  // Else set the damping coefficient:
  s.damping_coeff = 2. * oda_c <= -oda_s ? 1 : -oda_s / (2 * oda_c);
  assert_dbg(s.damping_coeff >= 0 && s.damping_coeff <= 1, krims::ExcInternalError());

#ifdef DEBUG_TODA
  // TODO temporary
  std::cout << std::setprecision(16) << "oda s:     " << oda_s << std::endl;
  std::cout << std::setprecision(16) << "oda c:     " << oda_c << std::endl;
  std::cout << std::setprecision(16) << "ST damping:  " << s.damping_coeff
            << std::setprecision(6) << std::endl;
#endif
}

template <typename ScfState>
void TruncatedOptDampScf<ScfState>::update_oda_diagmat(state_type& s) const {
  // Initialise an empty matrix for the tODA fock matrix guess in the state:
  s.diagonalised_matrix_ptr = std::make_shared<diagmat_type>();

  if (s.prev_problem_matrix_ptr == nullptr || s.damping_coeff == 1) {
    (*s.diagonalised_matrix_ptr) += s.problem_matrix();
  } else {
    (*s.diagonalised_matrix_ptr) += (1. - s.damping_coeff) * (*s.prev_problem_matrix_ptr);
    (*s.diagonalised_matrix_ptr) += s.damping_coeff * s.problem_matrix();
  }

  on_new_diagmat(s);
}

}  // namespace gscf
