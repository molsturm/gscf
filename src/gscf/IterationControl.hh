#pragma once
#include "IterationState.hh"
#include <cstddef>
#include <linalgwrap/Subscribable.hh>

namespace gscf {

/** Basic struct which contains parameters that control
 *  the way the iteration proceeds.
 *
 * \tparam IterationState The type of the iteration state
 *         to pass to is_converged.
 */
template <typename IterationState>
struct IterationControl : public linalgwrap::Subscribable {
  typedef IterationState it_state_type;
  typedef typename it_state_type::count_type count_type;

  static_assert(
        std::is_base_of<gscf::IterationState, IterationState>::value,
        "IterationTraits needs to be derived off gscf::IterationTraits");

  /** \brief Function to check convergence before a new iteration begins.
   *  This default implementation does no checking and just returns
   *  false.
   */
  virtual bool is_converged(const it_state_type&) const { return false; }

  /** Maximum number of iterations
   *
   * The iteration should be considered as failed
   * once we go beyond this number of iterations.
   **/
  count_type max_iter;

  /** Construct an IterationControl object
   *
   * \param max_iter_ Default maximum number of iterations.
   * */
  IterationControl(count_type max_iter_) : max_iter{max_iter_} {}
};

}  // namespace gscf
