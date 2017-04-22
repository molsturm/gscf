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

#include "pulay_error.hh"
#include "FockMatrix.hh"

namespace gscfmock {

matrix_type pulay_error(const FockMatrix& fock_bb,
                        const linalgwrap::MultiVector<const vector_type>& coefficients_bf,
                        const matrix_type& overlap_bb) {
  auto occ_a = fock_bb.indices_subspace(gscf::OrbitalSpace::OCC_ALPHA);
#ifdef DEBUG
  auto occ_b = fock_bb.indices_subspace(gscf::OrbitalSpace::OCC_BETA);
  assert_dbg(occ_a == occ_b, krims::ExcNotImplemented());
#endif

  // Occupied coefficients
  auto ca_bo = coefficients_bf.subview(occ_a);

  // Form first products (Factor of 2 since alpha == beta)
  auto Sca_bo = 2. * overlap_bb * ca_bo;
  auto Fca_bo = fock_bb * ca_bo;

  // Form the antisymmetric outer product and return it.
  // The idea is
  // S * P * F - F * P * S == S * C * C^T * F - F * C * C^T * S
  //                       == (S*C) * (F*C)^T - (F*C) * (S*C)^T
  return outer_prod_sum(Sca_bo, Fca_bo) - outer_prod_sum(Fca_bo, Sca_bo);
}

}  // namespace gscfmock
