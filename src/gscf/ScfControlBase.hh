#pragma once
#include "IterationControl.hh"
#include "ScfStateBase.hh"

namespace gscf {

/** \brief Class to provide the basic set of parameters and the convergence
 *  checks the ScfBase class needs.
 *
 *  Consider overriding the is_converged function to provide a convergence
 *  check.
 *
 * \tparam ScfState The type of the scf state. Via the type traits underlying
 *         the scf state the neccessary types for this class are also
 *         determined.
 *
 * \see gscf::ScfBase
 */
template <typename ScfState>
struct ScfControlBase : public IterationControl<ScfState> {
  typedef IterationControl<ScfState> base_type;
  typedef typename base_type::it_state_type scf_state_type;

  /** \name Types forwarded from ScfState */
  ///@{
  /** The type of the problem matrix as determined by the traits */
  typedef typename scf_state_type::probmat_type probmat_type;

  /** The type of the problem matrix as determined by the traits */
  typedef typename scf_state_type::diagmat_type diagmat_type;

  /** The type of the scalars as determined by the traits */
  typedef typename scf_state_type::scalar_type scalar_type;

  /** The type of the size indices as determined by the traits */
  typedef typename scf_state_type::size_type size_type;

  /** The type of the stored matrices as determined by the traits */
  typedef typename scf_state_type::matrix_type matrix_type;

  /** The type of the stored vectors as determined by the traits */
  typedef typename scf_state_type::vector_type vector_type;
  ///@}

  static_assert(IsScfState<ScfState>::value,
                "ScfState needs to be derived off ScfStateBase.");

  /** The number of eigenpairs to calculate in the SCF procedure */
  size_type n_eigenpairs;

  /** The key used to pass the new eigenvector coefficients to the
   * non-linear problem matrix within the update_problem_matrix step. */
  std::string update_key = "evec_coefficients";

  /** Construct an ScfControlBase object and set some defaults
   *
   * \param n_eigenpairs_    Number of eigenpairs to compute
   *                         in the SCF procedure
   * \param max_iter_        The maximum number of iterations to perform.
   * */
  ScfControlBase(size_type n_eigenpairs_, size_type max_iter_)
        : base_type(max_iter_), n_eigenpairs(n_eigenpairs_){};
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is derived from ScfControlBase
 *
 * The definition is done using SFINAE, such that even for types not having a
 * typedef scalar_type this expression is valid.
 *  */
template <typename T, typename = void>
struct IsScfControl : public std::false_type {};

template <typename T>
struct IsScfControl<T, linalgwrap::void_t<typename T::scf_state_type>>
      : public std::is_base_of<ScfControlBase<typename T::scf_state_type>, T> {
};
//@}

}  // gscf
