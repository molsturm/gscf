#pragma once
#include <gscf/PlainScf.hh>
#include <gscf/PulayDiisScf.hh>
#include <iostream>
#include <iterator>
#include <linalgwrap/io.hh>

namespace scf_dummy {
using namespace gscf;

template <typename Fock>
auto pulay_error(const Fock& fock_bb,
                 const typename Fock::stored_matrix_type& coefficients_bf,
                 const typename Fock::stored_matrix_type& overlap_bb) ->
      typename Fock::stored_matrix_type {
  typedef typename Fock::size_type size_type;
  typedef typename Fock::stored_matrix_type matrix_type;

  // Calculate Pulay scf error.
  using namespace linalgwrap;

  const size_type n_alpha = fock_bb.n_alpha();
#ifdef DEBUG
  const size_type n_beta = fock_bb.n_beta();
#endif
  assert_dbg(n_alpha == n_beta, linalgwrap::ExcNotImplemented());

  // Occupied coefficients
  auto ca_bo = view::columns(coefficients_bf, range(n_alpha));

  // Density: (Factor of 2 since alpha == beta)
  auto pa_bb = 2 * ca_bo * view::transpose(ca_bo);

  // Return error expression
  // == S * P * F - F * P * S
  return overlap_bb * static_cast<matrix_type>(pa_bb * fock_bb) -
         fock_bb * (pa_bb * overlap_bb);
}

template <typename ScfState>
struct PlainScfHartreeFockControl : public PlainScfControl<ScfState> {
  typedef PlainScfControl<ScfState> base_type;
  typedef typename base_type::scf_state_type scf_state_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  scalar_type max_frob_pulay_error;

  PlainScfHartreeFockControl() : base_type{}, max_frob_pulay_error{5e-7} {}

  /** Check convergence by checking the maxmial deviation of
   *  the last and previous eval pointers */
  bool is_converged(const scf_state_type& s) const override {
    if (s.problem_matrix_ptr() == nullptr || s.eigenvectors_ptr() == nullptr)
      return false;

    const probmat_type& fock_bb = *s.problem_matrix_ptr();
    const matrix_type& coefficients_bf = *s.eigenvectors_ptr();
    const matrix_type& overlap_bb = s.overlap_matrix();

    matrix_type error = pulay_error(fock_bb, coefficients_bf, overlap_bb);
    return error.norm_frobenius() < max_frob_pulay_error;
  }
};

template <typename FockType>
class PlainScfHartreeFock
      : public PlainScf<FockType, PlainScfState<FockType>,
                        PlainScfHartreeFockControl<PlainScfState<FockType>>> {
public:
  typedef FockType fock_type;
  typedef PlainScf<FockType, PlainScfState<FockType>,
                   PlainScfHartreeFockControl<PlainScfState<FockType>>>
        base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::scf_state_type scf_state_type;

  PlainScfHartreeFock(linalgwrap::io::DataWriter_i<scalar_type>& writer)
        : m_writer(writer) {}

protected:
  void before_iteration_step(scf_state_type& s) const override {
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
    auto error = pulay_error(problem_matrix, *s.eigenvectors_ptr(),
                             s.overlap_matrix());
    energies_type hf_energies = problem_matrix.energies();

    auto n_iter = s.n_iter_count();
    std::string itstr = std::to_string(n_iter);

    assert_dbg(problem_matrix.are_hf_terms_stored(),
               linalgwrap::ExcInvalidState(
                     "problem matrix does not store HF terms."));
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
              << "     pulay_error: " << error.norm_frobenius() << std::endl;
  }

private:
  linalgwrap::io::DataWriter_i<scalar_type>& m_writer;
};

template <typename FockType>
class DiisScfHartreeFock : public PulayDiisScf<FockType> {
public:
  typedef FockType fock_type;
  typedef PulayDiisScf<FockType> base_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;
  typedef typename base_type::scf_state_type scf_state_type;

  DiisScfHartreeFock(linalgwrap::io::DataWriter_i<scalar_type>& writer)
        : m_writer(writer) {}

protected:
  matrix_type calculate_error(const scf_state_type& s) const override {
    const fock_type& fock_bb = *s.problem_matrix_ptr();
    const matrix_type& coefficients_bf = *s.eigenvectors_ptr();
    const matrix_type& overlap_bb = s.overlap_matrix();
    return pulay_error(fock_bb, coefficients_bf, overlap_bb);
  }

  void before_iteration_step(scf_state_type& s) const override {
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

  void on_new_diis_diagmat(scf_state_type& s) const override {
    if (s.diis_coefficients.size() > 0) {
      std::cout << "   DIIS coefficients:   ";
      std::copy(s.diis_coefficients.begin(), s.diis_coefficients.end(),
                std::ostream_iterator<scalar_type>(std::cout, " "));
      std::cout << std::endl;
    }

    m_writer.write("diiscoeff" + std::to_string(s.n_iter_count()),
                   s.diis_coefficients);
    m_writer.write("diisdiagmat" + std::to_string(s.n_iter_count()),
                   *s.diagonalised_matrix_ptr());
  }

  void on_update_problem_matrix(scf_state_type& s) const override {
    typedef typename fock_type::energies_type energies_type;

    auto& problem_matrix = *s.problem_matrix_ptr();
    energies_type hf_energies = problem_matrix.energies();

    auto n_iter = s.n_iter_count();
    std::string itstr = std::to_string(n_iter);

    assert_dbg(problem_matrix.are_hf_terms_stored(),
               linalgwrap::ExcInvalidState(
                     "problem matrix does not store HF terms."));
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

  void after_iteration_step(scf_state_type& s) const override {
    std::cout << "   Current DIIS error: " << s.errors.back().norm_frobenius()
              << std::endl;
    m_writer.write("diiserror" + std::to_string(s.n_iter_count()),
                   s.errors.back());
    m_writer.write("diislinsysmat" + std::to_string(s.n_iter_count()),
                   base_type::diis_linear_system_matrix(s));
  }

private:
  linalgwrap::io::DataWriter_i<scalar_type>& m_writer;
};
}
