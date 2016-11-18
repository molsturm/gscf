#pragma once
#include <linalgwrap/Base/Solvers.hh>

namespace gscf {

/** \brief Scf state base class
 *
 * Provides the most basic state of each SCF, the current problem matrix,
 * the current eigenpairs.
 *
 * \tparam ProblemMatrix The type of the Problem matrix object.
 */
template <typename ProblemMatrix, typename DiagonalisedMatrix = ProblemMatrix>
class ScfStateBase : public linalgwrap::SolverStateBase,
                     public linalgwrap::IterativeSolverState {
public:
  /** \name Type definitions */
  ///@{
  /** The type of the problem matrix, i.e. the non-linear problem
   *  we which to solve.
   *
   *  We assume that the Problem is Hermitian (or symmetric).
   *  */
  typedef ProblemMatrix probmat_type;

  /** The type of the matrix to be diagonalised in order to obtain
   *  the new eigenpairs in the SCF algorithm. */
  typedef DiagonalisedMatrix diagmat_type;

  /** The type of the scalars */
  typedef typename probmat_type::scalar_type scalar_type;

  /** The type of the size indices as determined by the traits */
  typedef typename probmat_type::size_type size_type;

  /** The type of the stored matrices as determined by the traits */
  typedef typename probmat_type::stored_matrix_type matrix_type;

  /** The type of the eigenvectors */
  typedef typename matrix_type::vector_type vector_type;
  ///@}

  /** Get the overlap matrix of thc SCF problem */
  const matrix_type& overlap_matrix() const { return *m_overlap_matrix_ptr; }

  /** Constant access to the current problem matrix, i.e. the current
   * approximation
   *  to the self-consistent state we wish to obtain for the non-linear problem
   *  in the end. */
  const std::shared_ptr<const probmat_type> problem_matrix_ptr() const {
    return m_problem_matrix_ptr;
  }

  /** Access to the current problem matrix, i.e. the current approximation
   *  to the self-consistent state we wish to obtain for the non-linear problem
   *  in the end. */
  std::shared_ptr<probmat_type>& problem_matrix_ptr() {
    return m_problem_matrix_ptr;
  }

  /** \brief Constant access to the matrix which is diagonalised in order to
   * obtain
   * the new eigenpairs.
   *
   * \note This matrix may be identical to problem_matrix_ptr() for some
   * algorithms
   * (like PlainSCF), but may also me different (like in a DiisScf).
   */
  const std::shared_ptr<const diagmat_type> diagonalised_matrix_ptr() const {
    return m_diagonalised_matrix;
  }

  /** \brief Access to the matrix which is diagonalised in order to obtain
   * the new eigenpairs.
   *
   * \note This matrix may be identical to problem_matrix_ptr() for some
   * algorithms
   * (like PlainSCF), but may also me different (like in a DiisScf).
   */
  std::shared_ptr<diagmat_type>& diagonalised_matrix_ptr() {
    return m_diagonalised_matrix;
  }

  /** Constant access to the current eigenvectors */
  const std::shared_ptr<const linalgwrap::MultiVector<vector_type>>
  eigenvectors_ptr() const {
    return m_eigenvectors_ptr;
  }

  /** Access to the current eigenvectors */
  std::shared_ptr<linalgwrap::MultiVector<vector_type>>& eigenvectors_ptr() {
    return m_eigenvectors_ptr;
  }

  /** Constant access to the current eigenvalues */
  const std::shared_ptr<const std::vector<scalar_type>> eigenvalues_ptr()
        const {
    return m_eigenvalues_ptr;
  }

  /** Access to the current eigenvalues */
  std::shared_ptr<std::vector<scalar_type>>& eigenvalues_ptr() {
    return m_eigenvalues_ptr;
  }

  /** \brief Constructor
   *
   * Initialise the scf state from a problem matrix and a overlap matrix.
   *
   * All other entities are given some sensible results, e.g. the
   * eigenvectors and the eigenvalues pointers are set to nullptr.
   * */
  ScfStateBase(probmat_type prob_mat, const matrix_type& overlap_mat)
        : m_overlap_matrix_ptr{krims::make_subscription(overlap_mat,
                                                        "ScfState")},
          m_problem_matrix_ptr{
                std::make_shared<probmat_type>(std::move(prob_mat))},
          m_diagonalised_matrix{nullptr},
          m_eigenvectors_ptr{nullptr},
          m_eigenvalues_ptr{nullptr} {
    // Check that we really get a hermitian matrix as we implicitly assume.
    assert_dbg(m_problem_matrix_ptr->is_hermitian(),
               linalgwrap::ExcMatrixNotHermitian());
    // Note: This assumption is build into the update_eigenpairs method
    //       of ScfBase.
  }

private:
  //! The overlap matrix of the SCF problem.
  krims::SubscriptionPointer<const matrix_type> m_overlap_matrix_ptr;

  /** The current problem matrix as obtained as the approximation
   *  to the self-consistent problem matrix in this step. **/
  std::shared_ptr<probmat_type> m_problem_matrix_ptr;

  /** The current matrix used to obtain the eigenvectors
   *  and eigenpairs. (Note that this may be identical to
   *  m_problem_matrix_ptr like in a plain SCF, but may also be
   *  different like in a DIIS-SCF */
  std::shared_ptr<diagmat_type> m_diagonalised_matrix;

  //! The current set of eigenvectors of the operator
  std::shared_ptr<linalgwrap::MultiVector<vector_type>> m_eigenvectors_ptr;

  //! The current eigenvalues of the operator
  std::shared_ptr<std::vector<scalar_type>> m_eigenvalues_ptr;

  // Check that the eigenvalue and eigenvector types are as we expect.
  // We need to do this, since we want to keep the interface simple here.
  // The eigenpair method handles all kind of eigenproblems and returns
  // its findings in the most sensible type possible. For Hermitian problems
  // this is as of now always the vector_type of the stored matrix for
  // the eigenvectors and the scalar type for the eigenvalues.
  // This is just here to check that this is really what we get.
  static_assert(std::is_same<vector_type,
                             typename linalgwrap::EigensolutionTypeFor<
                                   true, diagmat_type>::evector_type>::value,
                "Problem when determining the eigenvector type");
  static_assert(std::is_same<scalar_type,
                             typename linalgwrap::EigensolutionTypeFor<
                                   true, diagmat_type>::evalue_type>::value,
                "Problem when determining the eigenvalue type");
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is derived from ScfStateBase
 *  */
template <typename T, typename = void>
struct IsScfState : public std::false_type {};

template <typename T>
struct IsScfState<
      T, krims::VoidType<typename T::probmat_type, typename T::diagmat_type>>
      : public std::is_base_of<
              ScfStateBase<typename T::probmat_type, typename T::diagmat_type>,
              T> {};
//@}

}  // gscf
