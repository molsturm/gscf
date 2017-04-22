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
#include <linalgwrap/MultiVector.hh>

namespace gscfmock {

template <typename Fock>
auto pulay_error(
      const Fock& fock_bb,
      const linalgwrap::MultiVector<const typename Fock::vector_type>& coefficients_bf,
      const typename Fock::stored_matrix_type& overlap_bb) ->
      typename Fock::stored_matrix_type {
  typedef typename Fock::size_type size_type;

  const size_type n_alpha = fock_bb.n_alpha();
#ifdef DEBUG
  const size_type n_beta = fock_bb.n_beta();
#endif
  assert_dbg(n_alpha == n_beta, krims::ExcNotImplemented());

  // Occupied coefficients
  auto ca_bo = coefficients_bf.subview(krims::range(n_alpha));

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
