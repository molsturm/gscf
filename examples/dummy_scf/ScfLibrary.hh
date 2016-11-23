#pragma once
#include <gscf/PlainScf.hh>
#include <gscf/PulayDiisScf.hh>
#include <gscfmock/pulay_error.hh>
#include <iostream>
#include <iterator>
#include <linalgwrap/io.hh>

namespace dummy_scf {
using namespace gscf;
using namespace linalgwrap;
using namespace gscfmock;

template <typename FockType>
class PlainScfHartreeFock
      : public PlainScf<
              PlainScfState<FockType, typename FockType::stored_matrix_type>> {
public:
  typedef FockType fock_type;
  typedef PlainScf<
        PlainScfState<FockType, typename FockType::stored_matrix_type>>
        base_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename matrix_type::vector_type vector_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::state_type state_type;

  //
  // Iteration control
  //
  scalar_type max_frob_pulay_error = 5e-7;

  /** Check convergence by checking the maxmial deviation of
   *  the last and previous eval pointers */
  bool is_converged(const state_type& s) const override {
    if (s.problem_matrix_ptr() == nullptr || s.eigenvectors_ptr() == nullptr)
      return false;

    const fock_type& fock_bb = *s.problem_matrix_ptr();
    const linalgwrap::MultiVector<vector_type>& coefficients_bf =
          *s.eigenvectors_ptr();
    const matrix_type& overlap_bb = s.overlap_matrix();

    matrix_type error = pulay_error(fock_bb, coefficients_bf, overlap_bb);
    return norm_frobenius(error) < max_frob_pulay_error;
  }

  //
  // Constructor
  // TODO not the way we would do it now
  // => we need a way to feed the solver with a ParameterMap
  PlainScfHartreeFock(linalgwrap::io::DataWriter_i<scalar_type>& writer)
        : m_writer(writer) {}

protected:
  void before_iteration_step(state_type& s) const override {
    std::cout << std::endl
              << "Starting SCF iteration " << s.n_iter_count() << std::endl;
  }

  void on_update_eigenpairs(state_type& s) const override {
    std::cout << "   New orbital eigenvalues: " << std::endl;

    assert_dbg(m_writer, krims::ExcIO());
    m_writer.write("evals" + std::to_string(s.n_iter_count()),
                   make_as_multivector<vector_type>(*s.eigenvalues_ptr()));
    m_writer.write("evecs" + std::to_string(s.n_iter_count()),
                   *s.eigenvectors_ptr());

    // Print orbital evals:
    std::ostream_iterator<scalar_type> out_it(std::cout, " ");
    std::cout << "        ";
    std::copy(s.eigenvalues_ptr()->begin(), s.eigenvalues_ptr()->end(), out_it);
    std::cout << std::endl;
  }

  void on_update_problem_matrix(state_type& s) const override {
    typedef typename fock_type::energies_type energies_type;

    auto& problem_matrix = *s.problem_matrix_ptr();
    auto error = pulay_error(problem_matrix, *s.eigenvectors_ptr(),
                             s.overlap_matrix());
    energies_type hf_energies = problem_matrix.energies();

    auto n_iter = s.n_iter_count();
    std::string itstr = std::to_string(n_iter);

    assert_dbg(
          problem_matrix.are_hf_terms_stored(),
          krims::ExcInvalidState("problem matrix does not store HF terms."));
    m_writer.write("pa" + itstr, problem_matrix.hf_terms().pa_bb);
    m_writer.write("pb" + itstr, problem_matrix.hf_terms().pb_bb);
    m_writer.write("j" + itstr, problem_matrix.hf_terms().j_bb);
    m_writer.write("ka" + itstr, problem_matrix.hf_terms().ka_bb);
    m_writer.write("ka" + itstr, problem_matrix.hf_terms().ka_bb);
    m_writer.write("kb" + itstr, problem_matrix.hf_terms().kb_bb);
    m_writer.write("fock" + itstr, problem_matrix);
    m_writer.write("error" + itstr, error);

    std::streamsize prec = std::cout.precision();
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
              << "     E_total  = " << std::setprecision(15)
              << hf_energies.energy_total << std::setprecision(prec)
              << std::endl
              << "     pulay_error: " << norm_frobenius(error) << std::endl;
  }

private:
  linalgwrap::io::DataWriter_i<scalar_type>& m_writer;
};

template <typename FockType>
class DiisScfHartreeFock
      : public PulayDiisScf<PulayDiisScfState<
              FockType, typename FockType::stored_matrix_type>> {
public:
  typedef FockType fock_type;
  typedef PulayDiisScf<
        PulayDiisScfState<FockType, typename FockType::stored_matrix_type>>
        base_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;
  typedef typename base_type::state_type state_type;

  DiisScfHartreeFock(linalgwrap::io::DataWriter_i<scalar_type>& writer)
        : m_writer(writer) {}

protected:
  matrix_type calculate_error(const state_type& s) const override {
    const fock_type& fock_bb = *s.problem_matrix_ptr();
    const MultiVector<vector_type>& coefficients_bf = *s.eigenvectors_ptr();
    const matrix_type& overlap_bb = s.overlap_matrix();
    return pulay_error(fock_bb, coefficients_bf, overlap_bb);
  }

  void before_iteration_step(state_type& s) const override {
    std::cout << std::endl
              << "Starting SCF iteration " << s.n_iter_count() << std::endl;
  }

  void on_update_eigenpairs(state_type& s) const override {
    std::cout << "   New orbital eigenvalues: " << std::endl;

    assert_dbg(m_writer, krims::ExcIO());
    m_writer.write("evals" + std::to_string(s.n_iter_count()),
                   make_as_multivector<vector_type>(*s.eigenvalues_ptr()));
    m_writer.write("evecs" + std::to_string(s.n_iter_count()),
                   *s.eigenvectors_ptr());

    // Print orbital evals:
    std::ostream_iterator<scalar_type> out_it(std::cout, " ");
    std::cout << "        ";
    std::copy(s.eigenvalues_ptr()->begin(), s.eigenvalues_ptr()->end(), out_it);
    std::cout << std::endl;
  }

  void on_new_diis_diagmat(state_type& s) const override {
    if (s.diis_coefficients.size() > 0) {
      std::cout << "   DIIS coefficients:   ";
      std::copy(s.diis_coefficients.begin(), s.diis_coefficients.end(),
                std::ostream_iterator<scalar_type>(std::cout, " "));
      std::cout << std::endl;
    }

    m_writer.write("diiscoeff" + std::to_string(s.n_iter_count()),
                   as_multivector(s.diis_coefficients));
    m_writer.write("diisdiagmat" + std::to_string(s.n_iter_count()),
                   *s.diagonalised_matrix_ptr());
  }

  void on_update_problem_matrix(state_type& s) const override {
    typedef typename fock_type::energies_type energies_type;

    auto& problem_matrix = *s.problem_matrix_ptr();
    energies_type hf_energies = problem_matrix.energies();

    auto n_iter = s.n_iter_count();
    std::string itstr = std::to_string(n_iter);

    assert_dbg(
          problem_matrix.are_hf_terms_stored(),
          krims::ExcInvalidState("problem matrix does not store HF terms."));
    m_writer.write("pa" + itstr, problem_matrix.hf_terms().pa_bb);
    m_writer.write("pb" + itstr, problem_matrix.hf_terms().pb_bb);
    m_writer.write("j" + itstr, problem_matrix.hf_terms().j_bb);
    m_writer.write("ka" + itstr, problem_matrix.hf_terms().ka_bb);
    m_writer.write("ka" + itstr, problem_matrix.hf_terms().ka_bb);
    m_writer.write("kb" + itstr, problem_matrix.hf_terms().kb_bb);
    m_writer.write("fock" + itstr, problem_matrix);

    std::streamsize prec = std::cout.precision();
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
              << "     E_total  = " << std::setprecision(15)
              << hf_energies.energy_total << std::setprecision(prec)
              << std::endl;
  }

  void after_iteration_step(state_type& s) const override {
    std::cout << "   Current DIIS error: " << norm_frobenius(s.errors.back())
              << std::endl;
    m_writer.write("diiserror" + std::to_string(s.n_iter_count()),
                   s.errors.back());
    m_writer.write("diislinsysmat" + std::to_string(s.n_iter_count()),
                   base_type::diis_linear_system_matrix(s));
  }

private:
  linalgwrap::io::DataWriter_i<scalar_type>& m_writer;
};
}  // namespace dummy_scf
