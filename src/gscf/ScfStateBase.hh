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
#include "ScfProblemMatrix_i.hh"
#include <lazyten/Base/Solvers.hh>

namespace gscf {

// TODO should this live in lazyten?
/** The type to access some eigenproblem statistics from the most recent
 *  inner eigensolver invocation. */
struct EigenproblemStatistics {
  /** Number of iterations the most recent eigensolver invocation needed
   *  to solve the problem. */
  size_t n_iter() const { return m_n_iter; }

  /** Number of matrix applies the most recent eigensolver invocation
   *  needed to solve the problem */
  size_t n_mtx_applies() const { return m_n_mtx_applies; }

  EigenproblemStatistics() {}
  EigenproblemStatistics(size_t n_iter, size_t n_mtx_applies)
        : m_n_iter(n_iter), m_n_mtx_applies(n_mtx_applies) {}

 private:
  size_t m_n_iter        = 0;
  size_t m_n_mtx_applies = 0;
};

/** \brief Scf state base class
 *
 * Provides the most basic state of each SCF, the current problem matrix,
 * the current eigenpairs.
 *
 * \tparam ProblemMatrix The type of the Problem matrix object.
 */
template <typename ProblemMatrix, typename OverlapMatrix,
          typename DiagonalisedMatrix = ProblemMatrix>
class ScfStateBase : public lazyten::IterativeStateWrapper<lazyten::SolverStateBase> {

  static_assert(std::is_base_of<ScfProblemMatrix_i, ProblemMatrix>::value,
                "The ProblemMatrix needs to be derived off ScfProblemMatrix_i for SCF "
                "algorithms to work.");

  static_assert(std::is_same<typename lazyten::StoredTypeOf<ProblemMatrix>::type,
                             typename lazyten::StoredTypeOf<OverlapMatrix>::type>::value,
                "The stored matrix types of OverlapMatrix and ProblemMatrix have to "
                "agree.");

 public:
  /** \name Type definitions */
  ///@{
  /** The type of the problem matrix, i.e. the non-linear problem
   *  we which to solve.
   *
   *  We assume that the Problem is Hermitian (or symmetric).
   *  */
  typedef ProblemMatrix probmat_type;

  /** The type of the overlap matrix */
  typedef OverlapMatrix overlap_type;

  /** The type of the matrix to be diagonalised in order to obtain
   *  the new eigenpairs in the SCF algorithm. */
  typedef DiagonalisedMatrix diagmat_type;

  /** The type of the scalars */
  typedef typename probmat_type::scalar_type scalar_type;

  /** The type for real values */
  typedef typename probmat_type::real_type real_type;

  /** The type of the size indices as determined by the traits */
  typedef typename probmat_type::size_type size_type;

  /** The type of the stored matrices as determined by the traits */
  typedef typename probmat_type::stored_matrix_type matrix_type;

  /** The type of the eigenvectors */
  typedef typename matrix_type::vector_type vector_type;

  /** The type of the eigensolution */
  typedef lazyten::Eigensolution<scalar_type, vector_type> esoln_type;
  ///@}

  /** Get the overlap matrix of the SCF problem */
  const OverlapMatrix& overlap_matrix() const { return *m_overlap_matrix_ptr; }

  /** Constant access to the current problem matrix, i.e. the current
   *  approximation to the self-consistent state we wish to obtain for
   *  the non-linear problem in the end. */
  const probmat_type& problem_matrix() const { return *problem_matrix_ptr; }

  /** Access to the current problem matrix, i.e. the current
   *  approximation to the self-consistent state we wish to obtain for
   *  the non-linear problem in the end. */
  probmat_type& problem_matrix() { return *problem_matrix_ptr; }

  /** \brief Constant access to the matrix which is diagonalised in order to
   * obtain the new eigenpairs.
   *
   * \note This matrix may be identical to problem_matrix_ptr() for some
   * algorithms (like PlainSCF), but may also me different (like in a
   * DiisScf).
   */
  const diagmat_type& diagonalised_matrix() const { return *diagonalised_matrix_ptr; }

  /** \brief Access to the matrix which is diagonalised in order to
   * obtain the new eigenpairs.
   *
   * \note This matrix may be identical to problem_matrix_ptr() for some
   * algorithms (like PlainSCF), but may also me different (like in a
   * DiisScf).
   */
  diagmat_type& diagonalised_matrix() { return *diagonalised_matrix_ptr; }

  /** Constant access to the most recent eigensolution obtained */
  const esoln_type& eigensolution() const { return m_eigensolution; }

  /** Access the eigensolution before the one stored in eigensolution() */
  const esoln_type& previous_eigensolution() const { return m_prev_eigensolution; }

  /** Access to some statistics (number of iterations, number of applies)
   *  of the most recent eigenproblem which was solved to obtain the
   *  eigensolution
   */
  const EigenproblemStatistics& eigenproblem_stats() const { return m_eprob_stats; }

  /** Update the eigensolution stored in the state along with its
   * statistics */
  void push_new_eigensolution(esoln_type new_eigensolution,
                              EigenproblemStatistics new_stats);

  /** The number of problem matrix applies needed so far to
   *  solve the scf problem
   */
  size_t n_mtx_applies() const override { return m_n_mtx_applies; }
  size_t& n_mtx_applies() { return m_n_mtx_applies; }

  /** Norm of the last error as computed by the SCF's calculate_error
   *  function */
  real_type last_error_norm;

  /** \brief Constructor
   *
   * Initialise the scf state from a problem matrix and a overlap matrix.
   *
   * All other entities are given some sensible results, e.g. the
   * eigenvectors and the eigenvalues pointers are set to nullptr.
   * */
  ScfStateBase(probmat_type prob_mat, const OverlapMatrix& overlap_mat)
        : lazyten::IterativeStateWrapper<
                lazyten::SolverStateBase>{lazyten::SolverStateBase()},
          last_error_norm{lazyten::Constants<real_type>::invalid},
          problem_matrix_ptr{std::make_shared<probmat_type>(std::move(prob_mat))},
          diagonalised_matrix_ptr{nullptr},
          m_n_mtx_applies(0),
          m_overlap_matrix_ptr{krims::make_subscription(overlap_mat, "ScfState")},
          m_eigensolution{},
          m_prev_eigensolution{} {
    // Check that we really get a hermitian matrix as we implicitly assume.
    // Note: This assumption is build into the update_eigenpairs method
    //       of ScfBase.
    assert_dbg(problem_matrix_ptr->is_hermitian(
                     100 * lazyten::Constants<real_type>::default_tolerance),
               lazyten::ExcMatrixNotHermitian());
  }

  /** \name Transfer a guess to this state */
  ///@{
  /** Setup the guess of this state by copying from another state.
   *
   * Also updates the state of the problem matrix to the guess solution
   * */
  template <typename DiagMat>
  void obtain_guess_from(const ScfStateBase<probmat_type, overlap_type, DiagMat>& other) {
    m_prev_eigensolution = other.previous_eigensolution();
    obtain_guess_from(other.eigensolution());
  }

  /** Setup the guess of this state by moving from another state.
   *
   * Also updates the state of the problem matrix to the guess solution
   * */
  template <typename DiagMat>
  void obtain_guess_from(ScfStateBase<probmat_type, overlap_type, DiagMat>&& other) {
    m_prev_eigensolution = std::move(other.m_prev_eigensolution);
    obtain_guess_from(std::move(other.m_eigensolution));
  }

  /** Setup the guess of this state by coping in another eigensolution.
   *
   * Also updates the state of the problem matrix to the guess solution
   */
  void obtain_guess_from(esoln_type other_soln);
  ///@}

  /** \name Advanced access to matrix pointers
   *
   * Use only if you know what you are doing.
   * */
  ///@{
  /** The current problem matrix as obtained as the approximation
   *  to the self-consistent problem matrix in this step. **/
  std::shared_ptr<probmat_type> problem_matrix_ptr;

  /** The current matrix used to obtain the eigenvectors
   *  and eigenpairs. (Note that this may be identical to
   *  m_problem_matrix_ptr like in a plain SCF, but may also be
   *  different like in a DIIS-SCF */
  std::shared_ptr<diagmat_type> diagonalised_matrix_ptr;
  ///@}
 private:
  /** Statistics about the most recent inner eigensolver call */
  EigenproblemStatistics m_eprob_stats;

  /** The number of operator applies in this solver */
  size_t m_n_mtx_applies;

  //! The overlap matrix of the SCF problem.
  krims::SubscriptionPointer<const overlap_type> m_overlap_matrix_ptr;

  //! The most recent eigensolution obtained (only contains pointers,
  //  so ok to store
  esoln_type m_eigensolution;

  //! The previous eigensolution, i.e. the one before m_eigensolution.
  esoln_type m_prev_eigensolution;

  // Check that the eigenvalue and eigenvector types are as we expect.
  // We need to do this, since we want to keep the interface simple here.
  // The eigenpair method handles all kind of eigenproblems and returns
  // its findings in the most sensible type possible. For Hermitian problems
  // this is as of now always the vector_type of the stored matrix for
  // the eigenvectors and the scalar type for the eigenvalues.
  // This is just here to check that this is really what we get.
  static_assert(
        std::is_same<esoln_type,
                     typename lazyten::EigensolutionTypeFor<true, diagmat_type>>::value,
        "Problem when determining the eigensolution type");
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is derived from ScfStateBase
 *  */
template <typename T, typename = void>
struct IsScfState : public std::false_type {};

template <typename T>
struct IsScfState<T, krims::VoidType<typename T::probmat_type, typename T::overlap_type,
                                     typename T::diagmat_type>>
      : public std::is_base_of<
              ScfStateBase<typename T::probmat_type, typename T::overlap_type,
                           typename T::diagmat_type>,
              T> {};
//@}

//
// --------------------------------------------------------
//

template <typename ProblemMatrix, typename OverlapMatrix, typename DiagonalisedMatrix>
void ScfStateBase<ProblemMatrix, OverlapMatrix, DiagonalisedMatrix>::
      push_new_eigensolution(esoln_type new_eigensolution,
                             EigenproblemStatistics new_stats) {
  // If memory is really an issue one should have the option to bin the previous
  // eigensolution instead.
  m_prev_eigensolution = std::move(m_eigensolution);
  m_eigensolution      = std::move(new_eigensolution);
  m_eprob_stats        = new_stats;
}

template <typename ProblemMatrix, typename OverlapMatrix, typename DiagonalisedMatrix>
void ScfStateBase<ProblemMatrix, OverlapMatrix, DiagonalisedMatrix>::obtain_guess_from(
      esoln_type other_soln) {
  m_eigensolution = std::move(other_soln);

  // Update problem matrix with this solution
  assert_dbg(
        m_eigensolution.evectors_ptr != nullptr &&
              m_eigensolution.evectors().n_vectors() > 0,
        krims::ExcInvalidState("Guess solution does not contain valid eigenvectors"));
  const std::string key = problem_matrix_ptr->scf_update_key();
  problem_matrix_ptr->update({{key, m_eigensolution.evectors_ptr}});
}

}  // gscf
