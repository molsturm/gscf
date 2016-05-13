#pragma once
#include "ScfBase.hh"

namespace gscf {

// TODO this stuff should go into a DIIS state
/** Number of previous steps we take into account for the DIIS.
 *
 * If ==
 */
// size_type n_prev_steps_considered;

//! Is DIIS accelleration on or off
// bool useDIIS;

//! Overlap between
// linalgwrap::SmallMatrix<scalar_type> error_overlaps;

// end TODO

template <typename ProblemMatrix>
class PulayDiisScfSolver : public ScfBase<ProblemMatrix> {
public:
  typedef ScfBase<ProblemMatrix> base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::scf_state_type scf_state_type;
  typedef typename base_type::vector_type vector_type;
};

//
// ------------------------------------------
//

template <typename ProblemMatrix>
typename PulayDiisScfSolver<ProblemMatrix>::scf_state_type
PulayDiisScfSolver<ProblemMatrix>::run(probmat_type probmat_bb,
                                       const matrix_type& metricmat_bb) const {
  size_t max_steps = 10;

  scf_state_type state(std::move(probmat_bb));
  for (size_t i = 0; i < max_steps; ++i) {
    base_type::start_iteration_step(state);
    base_type::solve_eigensystem(state, metricmat_bb);
    base_type::update_problem_matrix(state);

    if (state.is_failed() || state.is_converged()) break;
  }

  // Return the final state
  return state;
}
// calculate error --> abstract
// calculate error overlap --> implement
// caclulate diis coefficients
// form new guess

}  // gscf
