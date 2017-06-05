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
#include "Integrals/IntegralDataBase.hh"
#include <gscf/FocklikeMatrix_i.hh>
#include <linalgwrap/LazyMatrix_i.hh>

namespace gscfmock {

/** Struct to store the hf energies calculated in the FockMatrix class */
struct HFEnergies {
  //! The kinetic energy:
  scalar_type energy_kinetic;

  //! The electron - nuclear attraction energy:
  scalar_type energy_elec_nuc_attr;

  //! The energy of the one electron terms:
  scalar_type energy_1e_terms;

  //! The energy of the coulomb term:
  scalar_type energy_coulomb;

  //! The energy of the exchange term:
  scalar_type energy_exchange;

  //! The energy of the one electron terms:
  scalar_type energy_2e_terms;

  //! The energy of the one electron terms:
  scalar_type energy_total;

  HFEnergies()
        : energy_kinetic{0},
          energy_elec_nuc_attr(0),
          energy_1e_terms(0),
          energy_coulomb(0),
          energy_exchange(0),
          energy_2e_terms(0),
          energy_total(0) {}
};

struct HFTerms {
  //! The alpha density matrix for the current set of coefficients
  matrix_type pa_bb;

  //! The beta density matrix for the current set of coefficients
  matrix_type pb_bb;

  //! The coulomb matrix for the current set of coefficients
  matrix_type j_bb;

  //! The alpha exchange matrix for the current set of coefficients
  matrix_type ka_bb;

  //! The beta exchange matrix for the current set of coefficients
  matrix_type kb_bb;

  /** Initialise the components container
   *
   * \param nbas  Number of basis functions
   * \param zero  Initialise the values to zero.
   * */
  HFTerms(size_t nbas, bool zero = false)
        : pa_bb{nbas, nbas, zero},
          pb_bb{nbas, nbas, zero},
          j_bb{nbas, nbas, zero},
          ka_bb{nbas, nbas, zero},
          kb_bb{nbas, nbas, zero} {}
};

/** Very simple and slow implementation of a Fock matrix class for a
 * closed-shell System */
class FockMatrix final : public linalgwrap::LazyMatrix_i<matrix_type>,
                         public gscf::FocklikeMatrix_i<scalar_type> {
 public:
  typedef linalgwrap::LazyMatrix_i<matrix_type> base_type;
  typedef gscfmock::scalar_type scalar_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  /** Construct a Fock matrix from in integral_data object and an initial
   * guess, i.e. an initial set of coefficients.
   *
   * This class assumes that we calculate a ground state, closed shell
   * electron configuration for now.
   *
   * \param nalpha           Number of alpha electrons
   * \param nbeta            Number of beta electrons.
   * \param integral_data    The integral_data object
   * \param initial_guess_bf The initial guess matrix.
   *        If the provided basis set behind the IntegralData has
   *        nbas basis functions, then this matrix should have
   *        the dimensionality nbas x nfock, where nfock is the
   *        number of fock operator eigenstates we calculate.
   * \param store_hf_terms   Should the HF terms be stored
   *                         (exchange matrix, coulomb matrix, ...)
   */
  FockMatrix(size_t n_alpha, size_t n_beta, const IntegralDataBase& integral_data,
             const linalgwrap::MultiVector<vector_type>& initial_guess_bf,
             bool store_hf_terms = false);

  /** Return the number of rows of the matrix */
  size_t n_rows() const override { return m_fock_ptr->n_rows(); }

  /** Return the number of columns of the matrix */
  size_t n_cols() const override { return m_fock_ptr->n_cols(); }

  /** Return an element of the matrix */
  scalar_type operator()(size_t row, size_t col) const override {
    assert_greater(row, n_rows());
    assert_greater(col, n_cols());
    return (*m_fock_ptr)(row, col);
  }

  /** \brief Compute the Matrix-Multivector application
   *
   * Forward to stored fock matrix
   */
  template <typename VectorIn, typename VectorOut,
            linalgwrap::mat_vec_apply_enabled_t<FockMatrix, VectorIn, VectorOut>...>
  void apply(const linalgwrap::MultiVector<VectorIn>& x,
             linalgwrap::MultiVector<VectorOut>& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
             const scalar_type c_y = linalgwrap::Constants<scalar_type>::zero) const {
    m_fock_ptr->apply(x, y, mode, c_this, c_y);
  }

  /** \brief Compute the Matrix-Multivector application
   *
   * Forward to stored fock matrix
   */
  void apply(
        const linalgwrap::MultiVector<
              const linalgwrap::MutableMemoryVector_i<scalar_type>>& x,
        linalgwrap::MultiVector<linalgwrap::MutableMemoryVector_i<scalar_type>>& y,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_y = linalgwrap::Constants<scalar_type>::zero) const override {
    m_fock_ptr->apply(x, y, mode, c_this, c_y);
  }

  /** Clone the matrix */
  lazy_matrix_expression_ptr_type clone() const override {
    return lazy_matrix_expression_ptr_type(new FockMatrix(*this));
  }

  /** Return the scf update key */
  const std::string& scf_update_key() const override { return m_update_key; }

  /** Update the inner state:
   * Build the Fock matrix with the new coefficients
   *
   * It expects the new coefficients under the parameter key
   * returned by scf_update_key()
   */
  void update(const krims::GenMap& map) override {
    if (map.exists(m_update_key)) {
      build_fock_matrix_from_coefficient(
            map.at<const linalgwrap::MultiVector<vector_type>>(m_update_key));
    }
  }

  /** Return the current set of energies, calculated at the most recent update
   */
  const HFEnergies& energies() const { return m_energies; }

  scalar_type energy_1e_terms() const override { return m_energies.energy_1e_terms; }
  scalar_type energy_2e_terms() const override { return m_energies.energy_2e_terms; }

  /** Return the Hartree Fock terms
   *
   * asserts that they are actually stored in this class
   */
  const HFTerms& hf_terms() const {
    assert_dbg(m_terms_ptr, krims::ExcInvalidState("Terms not stored"));
    return *m_terms_ptr;
  }

  /** Are the Hartree Fock terms stored? */
  bool are_hf_terms_stored() const { return m_store_hf_terms; }

  virtual krims::Range<size_t> indices_orbspace(gscf::OrbitalSpace osp) const override;

 private:
  // Struct for the individual terms:

  /** Calculate the coulomb matrix from a given density */
  void calc_coulomb(const stored_matrix_type& density_bb,
                    stored_matrix_type& coul_bb) const;

  /** Calculate the exchange matrix from a given density */
  void calc_exchange(const stored_matrix_type& density_bb,
                     stored_matrix_type& exch_bb) const;

  /** Build the Fock matrix from a coefficient */
  void build_fock_matrix_from_coefficient(
        const linalgwrap::MultiVector<vector_type>& coefficients_bf);

  /** Build the Fock matrix from alpha and beta densities */
  void build_fock_matrix_from_density(stored_matrix_type pa_bb, stored_matrix_type pb_bb);

  /** Number of alpha electrons */
  size_t m_n_alpha;

  /** Number of beta electrons */
  size_t m_n_beta;

  //! The integral data
  krims::SubscriptionPointer<const IntegralDataBase> m_idata_ptr;

  //! The actual fock matrix for the current set of coefficients
  std::shared_ptr<stored_matrix_type> m_fock_ptr;

  //! Do we store the HF terms:
  bool m_store_hf_terms;

  /** Pointer to the SCF terms
   *
   * Only different from nullptr, when m_store_hf_terms is true.
   */
  std::shared_ptr<HFTerms> m_terms_ptr;

  //! The current energy values:
  HFEnergies m_energies;

  //! Key used for updating the state.
  const std::string m_update_key = "evec_coefficients";
};

}  // namespace gscfmock
