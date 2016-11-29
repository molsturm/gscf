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
template <typename ProblemMatrix, typename OverlapMatrix,
          typename DiagonalisedMatrix = ProblemMatrix>
class ScfStateBase
      : public linalgwrap::IterativeStateWrapper<linalgwrap::SolverStateBase> {

  static_assert(
        std::is_same<
              typename linalgwrap::StoredTypeOf<ProblemMatrix>::type,
              typename linalgwrap::StoredTypeOf<OverlapMatrix>::type>::value,
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
  typedef linalgwrap::Eigensolution<scalar_type, vector_type> esoln_type;
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
   * algorithms (like PlainSCF), but may also me different (like in a DiisScf).
   */
  const diagmat_type& diagonalised_matrix() const {
    return *diagonalised_matrix_ptr;
  }

  /** \brief Access to the matrix which is diagonalised in order to
   * obtain the new eigenpairs.
   *
   * \note This matrix may be identical to problem_matrix_ptr() for some
   * algorithms (like PlainSCF), but may also me different (like in a DiisScf).
   */
  diagmat_type& diagonalised_matrix() { return *diagonalised_matrix_ptr; }

  /** Constant access to the most recent eigensolution obtained */
  const esoln_type& eigensolution() const { return m_eigensolution; }

  /** Access to the most recent eigensolution obtained */
  esoln_type& eigensolution() { return m_eigensolution; }

  /** Number of iterations the most recent eigensolver invocation needed
   *  to solve the problem. */
  size_t n_eigenproblem_iter() const { return m_n_eigenproblem_iter; }
  size_t& n_eigenproblem_iter() { return m_n_eigenproblem_iter; }

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
        : linalgwrap::IterativeStateWrapper<
                linalgwrap::SolverStateBase>{linalgwrap::SolverStateBase()},
          last_error_norm{linalgwrap::Constants<real_type>::invalid},
          problem_matrix_ptr{
                std::make_shared<probmat_type>(std::move(prob_mat))},
          diagonalised_matrix_ptr{nullptr},
          m_n_eigenproblem_iter(0),
          m_overlap_matrix_ptr{
                krims::make_subscription(overlap_mat, "ScfState")},
          m_eigensolution{} {
    // Check that we really get a hermitian matrix as we implicitly assume.
    assert_dbg(problem_matrix_ptr->is_hermitian(),
               linalgwrap::ExcMatrixNotHermitian());
    // Note: This assumption is build into the update_eigenpairs method
    //       of ScfBase.
  }

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
  /** Number of iterations the most recent eigensolver took */
  size_t m_n_eigenproblem_iter;

  //! The overlap matrix of the SCF problem.
  krims::SubscriptionPointer<const overlap_type> m_overlap_matrix_ptr;

  //! The most recent eigensolution obtained (only contains pointers,
  //  so ok to store
  esoln_type m_eigensolution;

  // Check that the eigenvalue and eigenvector types are as we expect.
  // We need to do this, since we want to keep the interface simple here.
  // The eigenpair method handles all kind of eigenproblems and returns
  // its findings in the most sensible type possible. For Hermitian problems
  // this is as of now always the vector_type of the stored matrix for
  // the eigenvectors and the scalar type for the eigenvalues.
  // This is just here to check that this is really what we get.
  static_assert(
        std::is_same<esoln_type, typename linalgwrap::EigensolutionTypeFor<
                                       true, diagmat_type>>::value,
        "Problem when determining the eigensolution type");
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is derived from ScfStateBase
 *  */
template <typename T, typename = void>
struct IsScfState : public std::false_type {};

template <typename T>
struct IsScfState<
      T, krims::VoidType<typename T::probmat_type, typename T::overlap_type,
                         typename T::diagmat_type>>
      : public std::is_base_of<
              ScfStateBase<typename T::probmat_type, typename T::overlap_type,
                           typename T::diagmat_type>,
              T> {};
//@}

}  // gscf
