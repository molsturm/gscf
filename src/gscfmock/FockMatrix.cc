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

#include "FockMatrix.hh"
#include "Integrals.hh"

namespace gscfmock {

krims::Range<size_t> FockMatrix::indices_orbspace(gscf::OrbitalSpace osp) const {
  using gscf::OrbitalSpace;
  switch (osp) {
    case OrbitalSpace::OCC_ALPHA:
      return {0, m_n_alpha};
    case OrbitalSpace::OCC_BETA:
      return {0, m_n_alpha};
    case OrbitalSpace::VIRT_ALPHA:
      return {m_n_alpha, n_rows()};
    case OrbitalSpace::VIRT_BETA:
      return {m_n_beta, n_rows()};
  }

  assert_internal(false);
  return {0, 0};
}

FockMatrix::FockMatrix(size_t n_alpha, size_t n_beta,
                       const IntegralDataBase& integral_data,
                       const lazyten::MultiVector<vector_type>& initial_guess_bf,
                       bool store_hf_terms)
      : m_n_alpha{n_alpha},
        m_n_beta{n_beta},
        m_idata_ptr{krims::make_subscription(integral_data, "FockMatrix")},
        m_fock_ptr(std::make_shared<matrix_type>(integral_data.nbas(),
                                                 integral_data.nbas(), false)),
        m_store_hf_terms{store_hf_terms},
        m_terms_ptr{nullptr},
        m_energies() {
  build_fock_matrix_from_coefficient(initial_guess_bf);
}

void FockMatrix::calc_coulomb(const matrix_type& density_bb, matrix_type& coul_bb) const {
  assert_internal(density_bb.n_rows() == density_bb.n_cols());
  assert_internal(density_bb.n_rows() == coul_bb.n_rows());
  assert_internal(density_bb.n_cols() == coul_bb.n_cols());

  // Number of basis fuctions:
  size_t nbas = density_bb.n_rows();

  // J_{ab} = \sum_{cd} P_{cd} < ac | bd >
  // a and b are the same centre, so are c and d

  // Double loop over basis functions a and b:
  for (size_type a = 0; a < nbas; ++a) {
    for (size_type b = 0; b < nbas; ++b) {
      // Shell pair index for basis functions a and b:
      const size_type ab_pair = a * nbas + b;

      // Sum accumulator variable:
      scalar_type sum{0};

      // Double loop over basis functions c and d:
      for (size_type c = 0; c < nbas; ++c) {
        for (size_type d = 0; d < nbas; ++d) {
          // Shell pair index for basis functions c and d:
          const size_type cd_pair = c * nbas + d;

          // Perform contraction:
          sum += density_bb(c, d) * m_idata_ptr->i_bbbb()(ab_pair, cd_pair);
        }  // d
      }    // c

      coul_bb(a, b) = sum;
    }  // b
  }    // a
}

void FockMatrix::calc_exchange(const matrix_type& density_bb,
                               matrix_type& exch_bb) const {
  assert_internal(density_bb.n_rows() == density_bb.n_cols());

  // Number of basis fuctions:
  size_t nbas = density_bb.n_rows();

  // Zero the exchange matrix:
  exch_bb.set_zero();

  // K_{ab} = \sum_{cd} P_{cd} < ab | cd >
  // a and c are the same centre, so are b and d

  // Double loop over shell pair basis functions a and c
  for (size_type a = 0; a < nbas; ++a) {
    for (size_type c = 0; c < nbas; ++c) {
      // Shell pair index for basis functions a and c:
      const size_type ac_pair = a * nbas + c;

      // Double loop over shell pair basis functions d and b:
      for (size_type d = 0; d < nbas; ++d) {
        for (size_type b = 0; b < nbas; ++b) {
          // Shell pair index for d and b:
          const size_type db_pair = d * nbas + b;

          // Perform contraction:
          exch_bb(a, b) += density_bb(c, d) * m_idata_ptr->i_bbbb()(ac_pair, db_pair);
        }  // b
      }    // d
    }      // c
  }        // a
}

void FockMatrix::build_fock_matrix_from_density(matrix_type pa_bb, matrix_type pb_bb) {
  // Assert densities have the correct size:
  assert_size(m_idata_ptr->nbas(), pa_bb.n_rows());
  assert_size(m_idata_ptr->nbas(), pa_bb.n_cols());
  assert_size(m_idata_ptr->nbas(), pb_bb.n_rows());
  assert_size(m_idata_ptr->nbas(), pb_bb.n_cols());

  if (!m_terms_ptr || !m_terms_ptr.unique()) {
    // Either the term storage is empty or we are not the only
    // one using it, so allocate new storage:
    m_terms_ptr = std::make_shared<HFTerms>(m_idata_ptr->nbas(), false);
  }

  // Store densities:
  m_terms_ptr->pa_bb = std::move(pa_bb);
  m_terms_ptr->pb_bb = std::move(pb_bb);

  // Compute total density:
  matrix_type pt_bb = m_terms_ptr->pa_bb + m_terms_ptr->pb_bb;

  // Build 2e terms:
  calc_coulomb(pt_bb, m_terms_ptr->j_bb);
  calc_exchange(m_terms_ptr->pa_bb, m_terms_ptr->ka_bb);
  calc_exchange(m_terms_ptr->pb_bb, m_terms_ptr->kb_bb);

  // Get shortcut references
  const matrix_type& j_bb  = m_terms_ptr->j_bb;
  const matrix_type& ka_bb = m_terms_ptr->ka_bb;
  const matrix_type& kb_bb = m_terms_ptr->kb_bb;

  // Calculate 1e energy:
  m_energies.energy_kinetic       = trace(m_idata_ptr->t_bb() * pt_bb);
  m_energies.energy_elec_nuc_attr = trace(m_idata_ptr->v0_bb() * pt_bb);
  m_energies.energy_1e_terms =
        m_energies.energy_kinetic + m_energies.energy_elec_nuc_attr;

  // Calculate 2e energies and sum:
  m_energies.energy_coulomb  = 0.5 * trace(j_bb * pt_bb);
  m_energies.energy_exchange = -0.5 * trace(ka_bb * m_terms_ptr->pa_bb);
  m_energies.energy_exchange -= 0.5 * trace(kb_bb * m_terms_ptr->pb_bb);
  m_energies.energy_2e_terms = m_energies.energy_coulomb + m_energies.energy_exchange;
  m_energies.energy_total    = m_energies.energy_1e_terms + m_energies.energy_2e_terms;

  // If our fock pointer is not unique, i.e. there are other users,
  // make a new one first:
  if (!m_fock_ptr.unique()) {
    m_fock_ptr = std::make_shared<matrix_type>(m_fock_ptr->n_rows(), m_fock_ptr->n_cols(),
                                               false);
  }

  // Build two-electron terms and add everything to get fock matrix:
  assert_equal(m_n_alpha, m_n_beta);
  // TODO adapt this method and class for open-shell calculations!
  *m_fock_ptr = m_idata_ptr->t_bb()     // kinetic energy term
                + m_idata_ptr->v0_bb()  // nuc-electron attract. potential term
                + j_bb                  // coloumb of full density
                - ka_bb;                // exchange of same-spin density
  // TODO for beta fock matrix and uhf:
  //  m_fock_b = ... - build_exchange(pb_bb);

  // Deallocate storage if not needed any more:
  if (!m_store_hf_terms) {
    m_terms_ptr.reset();
  }
}

void FockMatrix::build_fock_matrix_from_coefficient(
      const lazyten::MultiVector<vector_type>& coefficients_bf) {
  // Assert coefficients have the correct size:
  assert_size(m_idata_ptr->nbas(), coefficients_bf.n_elem());

  // Views for the vectors of the occupied orbitals:
  auto ca_bo = coefficients_bf.subview(krims::range(m_n_alpha));
  auto cb_bo = coefficients_bf.subview(krims::range(m_n_beta));

  // Alpha and beta spin densities:
  auto pa_bb = outer_prod_sum(ca_bo, ca_bo);
  auto pb_bb = outer_prod_sum(cb_bo, cb_bo);

  build_fock_matrix_from_density(std::move(pa_bb), std::move(pb_bb));
}

}  // namespace gscfmock
