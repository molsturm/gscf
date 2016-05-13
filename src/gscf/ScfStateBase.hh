#pragma once
#include "IterationState.hh"
#include <linalgwrap/SubscriptionPointer.hh>
#include <linalgwrap/type_utils.hh>
#include <memory>

namespace gscf {

/** \brief Scf state base class
 *
 * Provides the most basic state of each SCF, the current problem matrix,
 * the current eigenpairs.
 */
template <typename ScfTraits>
class ScfStateBase : public IterationState {
public:
  typedef ScfTraits scf_traits_type;
  typedef typename scf_traits_type::probmat_type probmat_type;
  typedef typename scf_traits_type::scalar_type scalar_type;
  typedef typename scf_traits_type::size_type size_type;
  typedef typename scf_traits_type::matrix_type matrix_type;
  typedef typename scf_traits_type::vector_type vector_type;

  /** Get the metric matrix of thc SCF problem */
  const matrix_type metric_matrix() const { return *m_metric_matrix_ptr; }

  /** Constant access to the current problem matrix */
  const std::shared_ptr<const probmat_type> problem_matrix_ptr() const {
    return m_problem_matrix_ptr;
  }

  /** Access to the current problem matrix */
  std::shared_ptr<probmat_type>& problem_matrix_ptr() {
    return m_problem_matrix_ptr;
  }

  /** Constant access to the current eigenvector matrix */
  const std::shared_ptr<const matrix_type> eigenvectors_ptr() const {
    return m_eigenvectors_ptr;
  }

  /** Access to the current eigenvector matrix */
  std::shared_ptr<matrix_type>& eigenvectors_ptr() {
    return m_eigenvectors_ptr;
  }

  /** Constant access to the current eigenvalues */
  const std::shared_ptr<const vector_type> eigenvalues_ptr() const {
    return m_eigenvalues_ptr;
  }

  /** Access to the current eigenvalues */
  std::shared_ptr<vector_type>& eigenvalues_ptr() { return m_eigenvalues_ptr; }

  /** \brief Constructor
   *
   * Initialise the scf state from a problem matrix and a metric matrix.
   *
   * All other entities are given some sensible results, e.g. the
   * eigenvectors and the eigenvalues pointers are set to nullptr.
   * */
  ScfStateBase(probmat_type prob_mat, const matrix_type& metric_mat)
        : m_metric_matrix_ptr{linalgwrap::make_subscription(metric_mat,
                                                            "ScfState")},
          m_problem_matrix_ptr{
                std::make_shared<probmat_type>(std::move(prob_mat))},
          m_eigenvectors_ptr{nullptr},
          m_eigenvalues_ptr{nullptr} {}

private:
  //! The metric matrix of the SCF problem.
  linalgwrap::SubscriptionPointer<const matrix_type> m_metric_matrix_ptr;

  //! The current problem matrix
  std::shared_ptr<probmat_type> m_problem_matrix_ptr;

  //! The current set of eigenvectors of the operator
  std::shared_ptr<matrix_type> m_eigenvectors_ptr;

  //! The current eigenvalues of the operator
  std::shared_ptr<vector_type> m_eigenvalues_ptr;
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is derived from ScfStateBase
 *
 * The definition is done using SFINAE, such that even for types not having a
 * typedef scalar_type this expression is valid.
 *  */
template <typename T, typename = void>
struct IsScfState : public std::false_type {};

template <typename T>
struct IsScfState<T, linalgwrap::void_t<typename T::scf_traits_type>>
      : public std::is_base_of<ScfStateBase<typename T::scf_traits_type>, T> {};
//@}

}  // gscf
