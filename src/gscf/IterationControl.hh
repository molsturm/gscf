#pragma once
#include <cstddef>
#include <linalgwrap/Subscribable.hh>

namespace gscf {

/** Basic struct which contains parameters that control
 *  the way the iteration proceeds, including convergence
 *  parameters
 */
struct IterationControl : public linalgwrap::Subscribable {
  typedef size_t count_type;

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
