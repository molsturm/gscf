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

#include "loewdin_guess.hh"
#include <linalgwrap/eigensystem.hh>

namespace dummy_scf {
using namespace linalgwrap;

MultiVector<vector_type> loewdin_guess(const matrix_type& overlap_bb) {
  // apply Löwdin normalisation to the basis functions
  //   - Diagonalise the overlap
  //   - Take 1/\sqrt{evals} at the diagonal
  //   - results in orthonormalised basis functions

  try {
    auto sol = eigensystem_hermitian(overlap_bb);

    // Eigenvectors and eigenvalues.
    auto& evectors = sol.evectors();
    const auto& evalues = sol.evalues();

    assert_size(evectors.n_vectors(), overlap_bb.n_cols());
    assert_size(evectors.n_elem(), overlap_bb.n_rows());

    for (size_t i = 0; i < evectors.n_vectors(); ++i) {
      evectors[i] *= 1. / sqrt(evalues[i]);
    }

    return std::move(sol.evectors());
  } catch (const linalgwrap::SolverException& e) {
    std::cerr << "Obtaining Löwdin guess failed, using zero guess" << std::endl;
    return MultiVector<vector_type>(overlap_bb.n_rows(), overlap_bb.n_cols());
  }
}

}  // namespace dummy_scf
