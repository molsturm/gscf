#pragma once
#include "ScfBase.hh"
#include "ScfControlBase.hh"
#include "ScfStateBase.hh"
#include "ScfTraits.hh"
#include <linalgwrap/Constants.hh>
#include <linalgwrap/ParameterMap.hh>

namespace gscf {

template <typename ScfTraits>
struct PlainScfControl : public ScfControlBase<ScfTraits> {
    typedef ScfControlBase<ScfTraits> base_type;
    typedef typename base_type::probmat_type probmat_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::matrix_type matrix_type;
    typedef typename base_type::vector_type vector_type;

    //! Maximum difference in the coefficients
    scalar_type max_eigenvalue_change;

    PlainScfControl()
          : base_type{linalgwrap::Constants<size_type>::invalid, 100},
            max_eigenvalue_change{1e-6} {}
};

template <typename ScfTraits>
struct PlainScfState : public ScfStateBase<ScfTraits> {
    typedef ScfStateBase<ScfTraits> base_type;
    typedef typename base_type::probmat_type probmat_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::matrix_type matrix_type;
    typedef typename base_type::vector_type vector_type;

    std::shared_ptr<vector_type> previous_eigenvalues_ptr;

    PlainScfState(probmat_type probmat, const matrix_type& metric_mat)
          : base_type{std::move(probmat), metric_mat},
            previous_eigenvalues_ptr{nullptr} {}
};

template <typename ProblemMatrix>
class PlainScf : public ScfBase<PlainScfState<ScfTraits<ProblemMatrix>>> {
  public:
    typedef ScfBase<PlainScfState<ScfTraits<ProblemMatrix>>> base_type;
    typedef typename base_type::scf_state_type scf_state_type;
    typedef typename scf_state_type::scf_traits_type scf_traits_type;
    typedef PlainScfControl<scf_traits_type> scf_control_type;

    typedef typename scf_traits_type::probmat_type probmat_type;
    typedef typename scf_traits_type::scalar_type scalar_type;
    typedef typename scf_traits_type::size_type size_type;
    typedef typename scf_traits_type::matrix_type matrix_type;
    typedef typename scf_traits_type::vector_type vector_type;

    /** Run a Plain SCF eigenproblem with the provided problem matrix
     *  and metric matrix and return the final state.
     *
     * \note The returned state contains information on whether the SCF
     * converged or failed. If one wants to assert automatically that
     * the SCF converged, use solve_assert instead.
     *  */
    scf_state_type solve(probmat_type probmat_bb,
                         const matrix_type& metricmat_bb) const override;

    PlainScf(scf_control_type& scf_control);

  private:
    bool is_converged(const scf_state_type& s) const {
        // TODO this is a really bad convergence criterion ...

        if (!s.previous_eigenvalues_ptr) {
            // There are no prevous eigenvectors -> return false
            return false;
        }

        // Compute the difference:
        matrix_type difference =
              *s.eigenvalues_ptr() - *s.previous_eigenvalues_ptr;

        // Take absmax element of this diffenence
        // TODO Use matrix norm functions of linalgwrap once they are there
        scalar_type max = 0;
        for (auto elem : difference) {
            max = std::max(std::abs(elem), max);
        }

        if (max < m_scf_control_ptr->max_eigenvalue_change) return true;
        return false;
    }

    // TODO Notes about convergence to think about
    //    - We need a way to provide external convergence criteria, since
    //      we do not know about the HF energy
    //    - We need a mechanism to collect convergence information from
    //      multiple (internal and external) resources and overall
    //      decide based on this data whether to end or continue.
    //    - Maybe have an extra class that gets the state and just does
    //      convergence checking?
    //    - In any case the convergence check that is performed outside of
    //      this class should be conveyed via a virtual call and *not* a
    //      lambda (since the latter is manyfold slower)
    //
    //    It seems to me that this indicates that the basic iteration loop
    //    is probably better located in the base class and only a virtual
    //    iteration_step() function is really to overloaded and implemented
    //    in this class

    linalgwrap::SubscriptionPointer<scf_control_type> m_scf_control_ptr;
};

template <typename ProblemMatrix>
typename PlainScf<ProblemMatrix>::scf_state_type PlainScf<ProblemMatrix>::solve(
      probmat_type probmat_bb, const matrix_type& metricmat_bb) const {

    if (m_scf_control_ptr->n_eigenpairs ==
        linalgwrap::Constants<size_type>::invalid) {
        // Calculate them all:
        m_scf_control_ptr->n_eigenpairs = probmat_bb.n_rows();
    }

    // TODO Maybe do the loop inside ScfBase and just
    //      overload a pure virtual initialise, step, is_converged and
    //      finalise in here.

    scf_state_type state(std::move(probmat_bb), metricmat_bb);

    while (!state.is_failed() && !is_converged(state)) {
        base_type::start_iteration_step(state);

        // Update the pointer to the previous eigenvalues:
        state.previous_eigenvalues_ptr = state.eigenvalues_ptr();

        base_type::solve_eigensystem(state);
        base_type::update_problem_matrix(state);
        base_type::end_iteration_step(state);
    }

    // Return the final state
    return state;
}

template <typename ProblemMatrix>
PlainScf<ProblemMatrix>::PlainScf(scf_control_type& scf_control)
      : base_type{scf_control},
        m_scf_control_ptr{
              linalgwrap::make_subscription(scf_control, "PlainScf")} {}

}  // gscf
