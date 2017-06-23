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
#include "FocklikeMatrix_i.hh"
#include "OperatorLinearCombination.hh"
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
      : public ScfStateBase<ProblemMatrix, OverlapMatrix,
                            OperatorLinearCombination<ProblemMatrix>> {

  static_assert(IsFocklikeMatrix<ProblemMatrix>::value,
                "The ProblemMatrix has to be a fock-like matrix for this SCF algorithm.");

  typedef ScfStateBase<ProblemMatrix, OverlapMatrix,
                       OperatorLinearCombination<ProblemMatrix>>
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

  /** Fock prefactor threshold. If the prefactor in front
   *  of a Fock matrix becomes smaller than this value, than
   *  this Fock matrix is entirely ignored in order to reduce
   *  the number of matrix applies, which are performed.
   */
  real_type min_fock_prefactor = 0.05;

  /** Update control parameters from Parameter map */
  void update_control_params(const krims::GenMap& map) {
    base_type::update_control_params(map);
    n_prev_steps = map.at(TruncatedOptDampScfKeys::n_prev_steps, n_prev_steps);
    min_fock_prefactor =
          map.at(TruncatedOptDampScfKeys::min_fock_prefactor, min_fock_prefactor);
  }

  /** Get the current settings of all internal control parameters and
   *  update the GenMap accordingly.
   */
  void get_control_params(krims::GenMap& map) const {
    base_type::get_control_params(map);
    map.update(TruncatedOptDampScfKeys::n_prev_steps, n_prev_steps);
    map.update(TruncatedOptDampScfKeys::min_fock_prefactor, min_fock_prefactor);
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

 private:
  /** Compute the new damping coefficient from the values of s and c,
   * which we determined.
   */
  scalar_type compute_damping_coeff(scalar_type s, scalar_type c) const;

  /** Compute the trace $\tr F^{(n-1)} P^{(n)}$, where n is the current
   *  iteration count, i.e. s.problem_matrix() was built from the density
   *  $P^{(n))$ and s.prev_problem_matrix_ptr points to $F^{(n-1)}$.
   */
  scalar_type trace_fprev_dcur(state_type& s) const;
};

//
// ---------------------------------------------------------------
//

template <typename ScfState>
void TruncatedOptDampScf<ScfState>::solve_state(state_type& state) const {
  assert_dbg(!state.is_failed(), krims::ExcInvalidState("Cannot solve a failed state"));

  assert_implemented(n_prev_steps == 2);

  while (!base_type::convergence_reached(state)) {
    base_type::start_iteration_step(state);

    // Compute new damping coefficients
    update_damping_coefficients(state);
    update_oda_diagmat(state);

    // Solve the eigensystem of the tODA problem matrix
    base_type::update_eigenpairs(state);

    // Update the matrix applies, noting that there are
    // various terms in the diagmat matrix we actually diagonalise
    const size_t n_terms = state.damping_coeff == 1 ? 1ul : 2ul;
    state.n_mtx_applies() += n_terms * state.eigenproblem_stats().n_mtx_applies();

    // The ODA guess (stored in the diagonalised_matrix_ptr)
    // may contain copies of all operators in history.
    // To release the memory of what we don't need any more,
    // we reset it at this stage
    state.diagonalised_matrix_ptr.reset();

    {
      // Out of the two problem matrices we store keep the one which had the
      // largest contribution in the diagonalised_matrix. This thing is most
      // likely closer to the actual minimum.
      //
      // If we purge the more recent one, update the damping coefficient
      // accordingly.
      if (state.damping_coeff > 0.5) {
        state.prev_problem_matrix_ptr = state.problem_matrix_ptr;
      } else {
        state.damping_coeff = 1 - state.damping_coeff;
      }
    }

    // Form the next problem matrix from the eigensolution
    // determined above
    base_type::update_problem_matrix(state);

    // Calculate the new error.
    base_type::update_last_error_norm(state);

    base_type::end_iteration_step(state);
  }
}

template <typename ScfState>
typename TruncatedOptDampScf<ScfState>::scalar_type
TruncatedOptDampScf<ScfState>::trace_fprev_dcur(state_type& s) const {
  // The most recent eigensolution:
  const auto& coeff_bf = s.eigensolution().evectors();
  const auto& orben_f  = s.eigensolution().evalues();
  const auto& fp_bb    = *s.prev_problem_matrix_ptr;

  const krims::Range<size_t> occ_a = fp_bb.indices_orbspace(OrbitalSpace::OCC_ALPHA);
  const krims::Range<size_t> occ_b = fp_bb.indices_orbspace(OrbitalSpace::OCC_BETA);

  assert_greater_equal(occ_a.length(), orben_f.size());
  assert_greater_equal(occ_b.length(), orben_f.size());

  // In the special case where the damping coefficient is 1 the diagonalised
  // matrix equals the previous problem matrix exactly. Therefore we do not need to
  // perform the trace  and can just sum the eigenvalues of the occupied orbitals
  // to get the trace $\tr F^{(n)} P^{(n+1)}$
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
    //      block-diagonality in place.

    // increase apply count
    s.n_mtx_applies() += 1;

    // Occupied coefficients (alpha)
    const auto ca_bo = coeff_bf.subview(occ_a);

    // Apply fock to occupied (alpha) orbitals
    //     -- O(n_bas*n_bas*n_occ)
    const auto Fca_bo = fp_bb * ca_bo;

    // Form the outer product sum and trace it
    //  -- O(n_bas*n_bas*n_occ)
    const auto tr_aa = trace(outer_prod_sum(ca_bo, Fca_bo));

    if (occ_a == occ_b) {
      return 2. * tr_aa;
    } else {
      // Repeat the same in the beta blocks
      const auto cb_bo = coeff_bf.subview(occ_b);
      const auto tr_bb = trace(outer_prod_sum(cb_bo, fp_bb * cb_bo));
      return tr_aa + tr_bb;
    }
  }
}

template <typename ScfState>
void TruncatedOptDampScf<ScfState>::update_damping_coefficients(state_type& s) const {
  if (s.prev_problem_matrix_ptr == nullptr) {
    s.damping_coeff = 1;
    return;
  }

  // The value of the trace $\tr F^{(n-1)} P^{(n)}$
  const scalar_type tr_cpc = trace_fprev_dcur(s);

  // Compute s and c for the 2-step truncated ODA SCF
  const scalar_type oda_s = tr_cpc - s.prev_problem_matrix_ptr->energy_1e_terms() -
                            2. * s.prev_problem_matrix_ptr->energy_2e_terms();
  const scalar_type oda_c = s.problem_matrix().energy_1e_terms() - tr_cpc +
                            s.problem_matrix().energy_2e_terms() +
                            s.prev_problem_matrix_ptr->energy_2e_terms();

  s.damping_coeff = compute_damping_coeff(oda_s, oda_c);
  if (s.damping_coeff == 1) s.prev_problem_matrix_ptr.reset();
}

template <typename ScfState>
typename TruncatedOptDampScf<ScfState>::scalar_type
TruncatedOptDampScf<ScfState>::compute_damping_coeff(scalar_type oda_s,
                                                     scalar_type oda_c) const {
  // TODO Make threshold configurable
  // If the error in the SCF gets too small the numerics makes the values of
  // oda_s and oda_c go mad and it is not reliable to trust them at all.
  const real_type ignore_threshold = 1e-14;
  if (fabs(oda_s) < ignore_threshold || fabs(oda_c) < ignore_threshold) return 1;

  // Catch another frequent problem when numerical errors increase:
  if (oda_s >= 0) return 1;

  const scalar_type new_damping_coeff = 2. * oda_c <= -oda_s ? 1 : -oda_s / (2 * oda_c);
  assert_internal(new_damping_coeff >= 0 && new_damping_coeff <= 1);

  // Cap the value of the damping coefficient from below and above, such that
  // no fock prefactor gets smaller than min_fock_prefactor.
  // The rationale is that otherwise it's not worth doing two applies
  // in the next diagonalisation.
  if (new_damping_coeff < min_fock_prefactor) return min_fock_prefactor;
  if (new_damping_coeff > 1 - min_fock_prefactor) return 1;

  return new_damping_coeff;
}

template <typename ScfState>
void TruncatedOptDampScf<ScfState>::update_oda_diagmat(state_type& s) const {
  if (s.prev_problem_matrix_ptr == nullptr) {
    // Initialise with the problem matrix, since no other matrix known
    s.diagonalised_matrix_ptr = std::make_shared<diagmat_type>(s.problem_matrix());
  } else {
    // Form the linear combination between the new and the previous problem matrix
    s.diagonalised_matrix_ptr = std::make_shared<diagmat_type>(*s.prev_problem_matrix_ptr,
                                                               1. - s.damping_coeff);
    s.diagonalised_matrix_ptr->push_term(s.problem_matrix(), s.damping_coeff);
  }

  on_new_diagmat(s);
}

}  // namespace gscf
