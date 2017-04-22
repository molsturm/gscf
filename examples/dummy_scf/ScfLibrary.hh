//
// Copyright (C) 2017 by the gscf authors
//
// This file is part of gscf.
//
// gscf is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// gscf is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with gscf. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include <gscf/PlainScf.hh>
#include <gscf/PulayDiisScf.hh>
#include <gscf/TruncatedOptDampScf.hh>
#include <gscfmock/FockMatrix.hh>
#include <gscfmock/pulay_error.hh>
#include <iostream>
#include <iterator>
#include <linalgwrap/io.hh>

namespace dummy_scf {
using namespace gscf;
using namespace linalgwrap;
using namespace gscfmock;

class PlainScfHartreeFock : public PlainScf<PlainScfState<FockMatrix, matrix_type>> {
 public:
  //
  // Constructor
  // TODO not the way we would do it now
  // => we need a way to feed the solver with a GenMap
  PlainScfHartreeFock(linalgwrap::io::DataWriter_i<scalar_type>& writer)
        : m_writer(writer) {}

 protected:
  matrix_type calculate_error(const state_type& s) const override {
    return pulay_error(s.problem_matrix(), s.eigensolution().evectors(),
                       s.overlap_matrix());
  }

  void before_iteration_step(state_type& s) const override {
    std::cout << std::endl << "Starting SCF iteration " << s.n_iter() << std::endl;
  }

  void on_update_eigenpairs(state_type& s) const override {
    std::cout << "   New orbital eigenvalues: " << std::endl;

    assert_dbg(m_writer, krims::ExcIO());
    const auto& evalues = s.eigensolution().evalues();
    auto evectors = s.eigensolution().evectors();

    m_writer.write("evals" + std::to_string(s.n_iter()),
                   make_as_multivector<vector_type>(evalues));
    m_writer.write("evecs" + std::to_string(s.n_iter()), evectors);

    // Print orbital evals:
    std::ostream_iterator<scalar_type> out_it(std::cout, " ");
    std::cout << "        ";
    std::copy(evalues.begin(), evalues.end(), out_it);
    std::cout << std::endl;
  }

  void on_update_problem_matrix(state_type& s) const override {
    auto error = calculate_error(s);
    auto hf_energies = s.problem_matrix().energies();

    auto n_iter = s.n_iter();
    std::string itstr = std::to_string(n_iter);

    assert_dbg(s.problem_matrix().are_hf_terms_stored(),
               krims::ExcInvalidState("problem matrix does not store HF terms."));
    m_writer.write("pa" + itstr, s.problem_matrix().hf_terms().pa_bb);
    m_writer.write("pb" + itstr, s.problem_matrix().hf_terms().pb_bb);
    m_writer.write("j" + itstr, s.problem_matrix().hf_terms().j_bb);
    m_writer.write("ka" + itstr, s.problem_matrix().hf_terms().ka_bb);
    m_writer.write("ka" + itstr, s.problem_matrix().hf_terms().ka_bb);
    m_writer.write("kb" + itstr, s.problem_matrix().hf_terms().kb_bb);
    m_writer.write("fock" + itstr, s.problem_matrix());
    m_writer.write("error" + itstr, error);

    std::streamsize prec = std::cout.precision();
    std::cout << "   Current HF energies:" << std::endl
              << "     E_kin    = " << hf_energies.energy_kinetic << std::endl
              << "     E_v0     = " << hf_energies.energy_elec_nuc_attr << std::endl
              << "     E_coul   = " << hf_energies.energy_coulomb << std::endl
              << "     E_xchge  = " << hf_energies.energy_exchange << std::endl
              << std::endl
              << "     E_1e     = " << hf_energies.energy_1e_terms << std::endl
              << "     E_2e     = " << hf_energies.energy_2e_terms << std::endl
              << std::endl
              << "     E_total  = " << std::setprecision(15) << hf_energies.energy_total
              << std::setprecision(prec) << std::endl
              << "     pulay_error: " << norm_frobenius(error) << std::endl;
  }

 private:
  linalgwrap::io::DataWriter_i<scalar_type>& m_writer;
};

class DiisScfHartreeFock
      : public PulayDiisScf<PulayDiisScfState<FockMatrix, matrix_type>> {
 public:
  DiisScfHartreeFock(linalgwrap::io::DataWriter_i<scalar_type>& writer)
        : m_writer(writer) {}

 protected:
  matrix_type calculate_error(const state_type& s) const override {
    return pulay_error(s.problem_matrix(), s.eigensolution().evectors(),
                       s.overlap_matrix());
  }

  void before_iteration_step(state_type& s) const override {
    std::cout << std::endl << "Starting SCF iteration " << s.n_iter() << std::endl;
  }

  void on_update_eigenpairs(state_type& s) const override {
    std::cout << "   New orbital eigenvalues: " << std::endl;

    assert_dbg(m_writer, krims::ExcIO());
    const auto& evalues = s.eigensolution().evalues();
    auto evectors = s.eigensolution().evectors();

    m_writer.write("evals" + std::to_string(s.n_iter()),
                   make_as_multivector<vector_type>(evalues));
    m_writer.write("evecs" + std::to_string(s.n_iter()), evectors);

    // Print orbital evals:
    std::ostream_iterator<scalar_type> out_it(std::cout, " ");
    std::cout << "        ";
    std::copy(evalues.begin(), evalues.end(), out_it);
    std::cout << std::endl;
  }

  void on_new_diagmat(state_type& s) const override {
    if (s.diis_coefficients.size() > 0) {
      std::cout << "   DIIS coefficients:   ";
      std::copy(s.diis_coefficients.begin(), s.diis_coefficients.end(),
                std::ostream_iterator<scalar_type>(std::cout, " "));
      std::cout << std::endl;
    }

    m_writer.write("diiscoeff" + std::to_string(s.n_iter()),
                   as_multivector(s.diis_coefficients));
    m_writer.write("diisdiagmat" + std::to_string(s.n_iter()), s.diagonalised_matrix());
  }

  void on_update_problem_matrix(state_type& s) const override {
    auto hf_energies = s.problem_matrix().energies();
    auto n_iter = s.n_iter();
    std::string itstr = std::to_string(n_iter);

    assert_dbg(s.problem_matrix().are_hf_terms_stored(),
               krims::ExcInvalidState("problem matrix does not store HF terms."));
    m_writer.write("pa" + itstr, s.problem_matrix().hf_terms().pa_bb);
    m_writer.write("pb" + itstr, s.problem_matrix().hf_terms().pb_bb);
    m_writer.write("j" + itstr, s.problem_matrix().hf_terms().j_bb);
    m_writer.write("ka" + itstr, s.problem_matrix().hf_terms().ka_bb);
    m_writer.write("ka" + itstr, s.problem_matrix().hf_terms().ka_bb);
    m_writer.write("kb" + itstr, s.problem_matrix().hf_terms().kb_bb);
    m_writer.write("fock" + itstr, s.problem_matrix());

    std::streamsize prec = std::cout.precision();
    std::cout << "   Current HF energies:" << std::endl
              << "     E_kin    = " << hf_energies.energy_kinetic << std::endl
              << "     E_v0     = " << hf_energies.energy_elec_nuc_attr << std::endl
              << "     E_coul   = " << hf_energies.energy_coulomb << std::endl
              << "     E_xchge  = " << hf_energies.energy_exchange << std::endl
              << std::endl
              << "     E_1e     = " << hf_energies.energy_1e_terms << std::endl
              << "     E_2e     = " << hf_energies.energy_2e_terms << std::endl
              << std::endl
              << "     E_total  = " << std::setprecision(15) << hf_energies.energy_total
              << std::setprecision(prec) << std::endl;
  }

  void after_iteration_step(state_type& s) const override {
    std::cout << "   Current DIIS error: " << norm_frobenius(s.errors.back())
              << std::endl;
    m_writer.write("diiserror" + std::to_string(s.n_iter()), s.errors.back());
    m_writer.write("diislinsysmat" + std::to_string(s.n_iter()),
                   diis_linear_system_matrix(s));
  }

 private:
  linalgwrap::io::DataWriter_i<scalar_type>& m_writer;
};

class TODAScfHartreeFock
      : public TruncatedOptDampScf<TruncatedOptDampScfState<FockMatrix, matrix_type>> {
 public:
  TODAScfHartreeFock(linalgwrap::io::DataWriter_i<scalar_type>& writer)
        : m_writer(writer) {}

 protected:
  matrix_type calculate_error(const state_type& s) const override {
    return pulay_error(s.problem_matrix(), s.eigensolution().evectors(),
                       s.overlap_matrix());
  }

  void before_iteration_step(state_type& s) const override {
    std::cout << std::endl << "Starting SCF iteration " << s.n_iter() << std::endl;
  }

  void on_update_eigenpairs(state_type& s) const override {
    std::cout << "   New orbital eigenvalues: " << std::endl;

    assert_dbg(m_writer, krims::ExcIO());
    const auto& evalues = s.eigensolution().evalues();
    auto evectors = s.eigensolution().evectors();

    m_writer.write("evals" + std::to_string(s.n_iter()),
                   make_as_multivector<vector_type>(evalues));
    m_writer.write("evecs" + std::to_string(s.n_iter()), evectors);

    // Print orbital evals:
    std::ostream_iterator<scalar_type> out_it(std::cout, " ");
    std::cout << "        ";
    std::copy(evalues.begin(), evalues.end(), out_it);
    std::cout << std::endl;
  }

  void on_new_diagmat(state_type& s) const override {
    if (s.prev_problem_matrix_ptr != nullptr) {
      std::cout << "   ODA damping coefficients:   ";
      std::cout << s.damping_coeff << std::endl;
    }

    m_writer.write("todacoeff" + std::to_string(s.n_iter()), s.damping_coeff);
    m_writer.write("todadiagmat" + std::to_string(s.n_iter()), s.diagonalised_matrix());
  }

  void on_update_problem_matrix(state_type& s) const override {
    auto hf_energies = s.problem_matrix().energies();
    auto n_iter = s.n_iter();
    std::string itstr = std::to_string(n_iter);

    assert_dbg(s.problem_matrix().are_hf_terms_stored(),
               krims::ExcInvalidState("problem matrix does not store HF terms."));
    m_writer.write("pa" + itstr, s.problem_matrix().hf_terms().pa_bb);
    m_writer.write("pb" + itstr, s.problem_matrix().hf_terms().pb_bb);
    m_writer.write("j" + itstr, s.problem_matrix().hf_terms().j_bb);
    m_writer.write("ka" + itstr, s.problem_matrix().hf_terms().ka_bb);
    m_writer.write("ka" + itstr, s.problem_matrix().hf_terms().ka_bb);
    m_writer.write("kb" + itstr, s.problem_matrix().hf_terms().kb_bb);
    m_writer.write("fock" + itstr, s.problem_matrix());

    std::streamsize prec = std::cout.precision();
    std::cout << "   Current HF energies:" << std::endl
              << "     E_kin    = " << hf_energies.energy_kinetic << std::endl
              << "     E_v0     = " << hf_energies.energy_elec_nuc_attr << std::endl
              << "     E_coul   = " << hf_energies.energy_coulomb << std::endl
              << "     E_xchge  = " << hf_energies.energy_exchange << std::endl
              << std::endl
              << "     E_1e     = " << hf_energies.energy_1e_terms << std::endl
              << "     E_2e     = " << hf_energies.energy_2e_terms << std::endl
              << std::endl
              << "     E_total  = " << std::setprecision(15) << hf_energies.energy_total
              << std::setprecision(prec) << std::endl;
  }

  void after_iteration_step(state_type& s) const override {
    std::cout << "   Current Pulay error: " << s.last_error_norm << std::endl;

    std::string itstr = std::to_string(s.n_iter());
    m_writer.write("pulayerror" + itstr, s.last_error_norm);
  }

 private:
  linalgwrap::io::DataWriter_i<scalar_type>& m_writer;
};

}  // namespace dummy_scf
