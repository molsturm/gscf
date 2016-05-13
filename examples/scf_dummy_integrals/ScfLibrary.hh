#pragma once
#include <gscf/PlainScf.hh>
#include <iostream>
#include <iterator>
#include <linalgwrap/io.hh>

namespace scf_dummy {
using namespace gscf;

template <typename ScfTraits>
struct PlainScfHartreeFockState : public PlainScfState<ScfTraits> {
  typedef PlainScfState<ScfTraits> base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  std::shared_ptr<vector_type> last_step_eval_ptr;

  PlainScfHartreeFockState(probmat_type probmat, const matrix_type& overlap_mat)
        : base_type{std::move(probmat), overlap_mat},
          last_step_eval_ptr{nullptr} {}
};

template <typename ScfState>
struct PlainScfHartreeFockControl : public PlainScfControl<ScfState> {
  typedef PlainScfControl<ScfState> base_type;
  typedef typename base_type::scf_state_type scf_state_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  scalar_type max_eval_difference;

  PlainScfHartreeFockControl() : base_type{}, max_eval_difference{1e-6} {}

  /** Check convergence by checking the maxmial deviation of
   *  the last and previous eval pointers */
  bool is_converged(const scf_state_type& state) const override {
    // If no previous eval_ptr, then we are not converged:
    if (!state.last_step_eval_ptr) return false;

    // Difference between eigenvalues:
    vector_type diff = *state.eigenvalues_ptr() - *state.last_step_eval_ptr;

    // Check difference is not too large for convergence:
    for (auto diff_elem : diff) {
      if (std::fabs(diff_elem) > max_eval_difference) return false;
    }
    return true;
  }
};

template <typename FockType>
class PlainScfHartreeFock
      : public PlainScf<FockType, PlainScfHartreeFockState<FockType>> {
public:
  typedef FockType fock_type;
  typedef PlainScf<FockType, PlainScfHartreeFockState<FockType>> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::scf_state_type scf_state_type;
  typedef PlainScfHartreeFockControl<scf_state_type> scf_control_type;

  PlainScfHartreeFock(scf_control_type& scf_control,
                      linalgwrap::io::DataWriter_i<scalar_type>& writer)
        : base_type(scf_control), m_writer(writer) {}

protected:
  void before_iteration_step(scf_state_type& s) const override {
    // Store the previous eigenvalues away for use in the convergence check:
    s.last_step_eval_ptr = s.eigenvalues_ptr();

    std::cout << std::endl
              << "Starting SCF iteration " << s.n_iter_count() << std::endl;
  }

  void on_update_eigenpairs(scf_state_type& s) const override {
    std::cout << "   New orbital eigenvalues: " << std::endl;

    assert_dbg(m_writer, linalgwrap::ExcIO());
    m_writer.write("evals" + std::to_string(s.n_iter_count()),
                   *s.eigenvalues_ptr());
    m_writer.write("evecs" + std::to_string(s.n_iter_count()),
                   *s.eigenvectors_ptr());

    // Print orbital evals:
    std::ostream_iterator<scalar_type> out_it(std::cout, " ");
    std::cout << "        ";
    std::copy(s.eigenvalues_ptr()->begin(), s.eigenvalues_ptr()->end(), out_it);
    std::cout << std::endl;
  }

  void on_update_problem_matrix(scf_state_type& s) const override {
    typedef typename fock_type::energies_type energies_type;

    auto& problem_matrix = *s.problem_matrix_ptr();
    energies_type hf_energies = problem_matrix.energies();

    auto n_iter = s.n_iter_count();
    std::string itstr = std::to_string(n_iter);

    assert_dbg(problem_matrix.are_hf_terms_stored(),
               ExcInvalidState("problem matrix does not store HF terms."));
    m_writer.write("pa" + itstr, problem_matrix.hf_terms().pa_bb);
    m_writer.write("pb" + itstr, problem_matrix.hf_terms().pb_bb);
    m_writer.write("j" + itstr, problem_matrix.hf_terms().j_bb);
    m_writer.write("ka" + itstr, problem_matrix.hf_terms().ka_bb);
    m_writer.write("ka" + itstr, problem_matrix.hf_terms().ka_bb);
    m_writer.write("kb" + itstr, problem_matrix.hf_terms().kb_bb);
    m_writer.write("fock" + itstr, problem_matrix);

    std::cout << "   Current HF energies:" << std::endl
              << "     E_kin    = " << hf_energies.energy_kinetic << std::endl
              << "     E_v0     = " << hf_energies.energy_elec_nuc_attr
              << std::endl
              << "     E_coul   = " << hf_energies.energy_coulomb << std::endl
              << "     E_xchge  = " << hf_energies.energy_exchange << std::endl
              << std::endl
              << "     E_1e     = " << hf_energies.energy_1e_terms << std::endl
              << "     E_2e     = " << hf_energies.energy_2e_terms << std::endl
              << std::endl
              << "     E_total  = " << hf_energies.energy_total << std::endl;
  }

private:
  linalgwrap::io::DataWriter_i<scalar_type>& m_writer;
};
}
