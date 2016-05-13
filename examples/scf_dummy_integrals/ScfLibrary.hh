#pragma once
#include <gscf/PlainScf.hh>
#include <iostream>
#include <iterator>
#include <linalgwrap/io.hh>

namespace scf_dummy {
using namespace gscf;

template <typename FockType>
class PlainScfHartreeFock : public PlainScf<FockType> {
public:
  typedef FockType fock_type;
  typedef PlainScf<FockType> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::scf_state_type scf_state_type;
  typedef typename base_type::scf_control_type scf_control_type;

  PlainScfHartreeFock(scf_control_type& scf_control,
                      linalgwrap::io::DataWriter_i<scalar_type>& writer)
        : base_type(scf_control), m_writer(writer) {}

protected:
  void before_iteration_step(const scf_state_type& s) const override {
    std::cout << std::endl
              << "Starting SCF iteration " << s.n_iter_count() << std::endl;
  }

  void on_update_eigenpairs(const scf_state_type& s) const override {
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

  void on_update_problem_matrix(const scf_state_type& s) const override {
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
