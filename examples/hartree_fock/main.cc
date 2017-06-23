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

// Setup the krims exception system for the tests.
#define KRIMS_INIT_EXCEPTION_SYSTEM
#include <krims/ExceptionSystem.hh>

#include "ScfLibrary.hh"
#include "loewdin_guess.hh"
#include <gscf/version.hh>
#include <gscfmock/FockMatrix.hh>
#include <gscfmock/Integrals.hh>
#include <linalgwrap/version.hh>

namespace dummy_scf {
/** Run an (atomic) SCF based on the integral data Sturmian14.
 *
 * \param nuc_charge   Nuclear charge
 * \param k_exp        Sturmian exponent
 * \param n_alpha      Number of alpha electrons
 * \param n_beta       Number of beta electrons
 */
void run_sturmian14(double nuc_charge, double k_exp, size_t n_alpha, size_t n_beta,
                    const std::string& method) {
  using gscfmock::IntegralsSturmian14;
  using gscfmock::FockMatrix;

  std::cout << "##############################################\n"
            << "#--  " << method << "  --#\n"
            << "#####" << std::string(method.size(), '#') << "#####" << std::endl;

  // Output file for mathematica debug:
  std::stringstream filename;
  filename << "/tmp/debug_gscf_scfdummy_" << method << "_" << nuc_charge << "_" << k_exp
           << "_sturm14.m";
  std::ofstream mathematicafile(filename.str());
  auto debugout = linalgwrap::io::make_writer<linalgwrap::io::Mathematica, scalar_type>(
        mathematicafile, 1e-5);

  // Define integral data and obtain a guess
  IntegralsSturmian14 idata(nuc_charge, k_exp);
  auto guess = loewdin_guess(idata.s_bb());

  // Write Sbb and guess to debug:
  debugout.write("sbb", idata.s_bb());
  debugout.write("guess", guess);

  // Setup Fock
  bool store_terms = true;
  FockMatrix fock{n_alpha, n_beta, idata, guess, store_terms};
  debugout.write("guessfock", fock);

  // Solve SCF
  if (method == "plain") {
    PlainScfHartreeFock{debugout}.solve(fock, idata.s_bb());
  } else if (method == "diis") {
    DiisScfHartreeFock{debugout}.solve(fock, idata.s_bb());
  } else if (method == "toda") {
    TODAScfHartreeFock{debugout}.solve(fock, idata.s_bb());
  }

  std::cout << "##############################################" << std::endl;
}
}  // namespace dummy_scf

int main() {
  const std::string atom  = "Be";
  const double k_exp      = 1.351;
  const size_t nuc_charge = 4.;
  const size_t n_alpha    = nuc_charge / 2;
  const size_t n_beta     = nuc_charge / 2;

  std::cout << "gscf version: " << gscf::version::version_string() << std::endl
            << "linalgwrap version: " << linalgwrap::version::version_string()
            << std::endl
            << "We compute closed-shell " << atom << " atom using k_exp == " << k_exp
            << std::endl
            << std::endl;

  for (const auto& method : {"plain", "diis", "toda"}) {
    dummy_scf::run_sturmian14(nuc_charge, k_exp, n_alpha, n_beta, method);
  }

  return 0;
}
