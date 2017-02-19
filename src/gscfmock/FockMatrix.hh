#pragma once
#include "Integrals/IntegralDataBase.hh"
#include <linalgwrap/LazyMatrix_i.hh>
#include <linalgwrap/io.hh>

namespace gscfmock {

/** Struct to store the hf energies calculated in the FockMatrix class */
template <typename Scalar>
struct HFEnergies {
  typedef Scalar scalar_type;

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
};

template <typename StoredMatrix>
struct HFTerms {
  typedef StoredMatrix stored_matrix_type;
  typedef typename stored_matrix_type::size_type size_type;

  //! The alpha density matrix for the current set of coefficients
  stored_matrix_type pa_bb;

  //! The beta density matrix for the current set of coefficients
  stored_matrix_type pb_bb;

  //! The coulomb matrix for the current set of coefficients
  stored_matrix_type j_bb;

  //! The alpha exchange matrix for the current set of coefficients
  stored_matrix_type ka_bb;

  //! The beta exchange matrix for the current set of coefficients
  stored_matrix_type kb_bb;

  /** Initialise the components container
   *
   * \param nbas  Number of basis functions
   * \param zero  Initialise the values to zero.
   * */
  HFTerms(size_type nbas, bool zero = false)
        : pa_bb{nbas, nbas, zero},
          pb_bb{nbas, nbas, zero},
          j_bb{nbas, nbas, zero},
          ka_bb{nbas, nbas, zero},
          kb_bb{nbas, nbas, zero} {}
};

/** Very simple and slow implementation of a Fock matrix class for a
 * closed-shell System */
template <typename IntegralData>
class FockMatrix : public linalgwrap::LazyMatrix_i<typename IntegralData::matrix_type> {
  static_assert(std::is_base_of<IntegralDataBase, IntegralData>::value,
                "IntegralData needs to be derived from IntegralDataBase");

 public:
  typedef IntegralData idata_type;
  typedef linalgwrap::LazyMatrix_i<typename idata_type::matrix_type> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::vector_type vector_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;
  typedef HFEnergies<scalar_type> energies_type;
  typedef HFTerms<stored_matrix_type> terms_type;

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
  FockMatrix(size_type n_alpha, size_type n_beta, const idata_type& integral_data,
             const linalgwrap::MultiVector<vector_type>& initial_guess_bf,
             bool store_hf_terms = false);

  /** Return the number of rows of the matrix */
  size_type n_rows() const override { return m_fock_ptr->n_rows(); }

  /** Return the number of columns of the matrix */
  size_type n_cols() const override { return m_fock_ptr->n_cols(); }

  /** Return an element of the matrix */
  scalar_type operator()(size_type row, size_type col) const override {
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
    return lazy_matrix_expression_ptr_type(new FockMatrix<IntegralData>(*this));
  }

  /** Return the scf update key */
  const std::string& scf_update_key() const { return m_update_key; }

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
  const energies_type& energies() const { return m_energies; }

  /** Return the Hartree Fock terms
   *
   * asserts that they are actually stored in this class
   */
  const terms_type& hf_terms() const {
    assert_dbg(m_terms_ptr, krims::ExcInvalidState("Terms not stored"));
    return *m_terms_ptr;
  }

  /** Are the Hartree Fock terms stored? */
  bool are_hf_terms_stored() const { return m_store_hf_terms; }

  /** Get the number of alpha electrons */
  size_type n_alpha() const { return m_n_alpha; }

  // TODO being able to access the above and below here feels wrong.
  // but we need it for the error calculation.

  /** Get the number of beta electrons */
  size_type n_beta() const { return m_n_beta; }

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
  size_type m_n_alpha;

  /** Number of beta electrons */
  size_type m_n_beta;

  //! The integral data
  krims::SubscriptionPointer<const idata_type> m_idata_ptr;

  //! The actual fock matrix for the current set of coefficients
  std::shared_ptr<stored_matrix_type> m_fock_ptr;

  //! Do we store the HF terms:
  bool m_store_hf_terms;

  /** Pointer to the SCF terms
   *
   * Only different from nullptr, when m_store_hf_terms is true.
   */
  std::shared_ptr<terms_type> m_terms_ptr;

  //! The current energy values:
  energies_type m_energies;

  //! Key used for updating the state.
  const std::string m_update_key = "evec_coefficients";
};

//
// --------------------------------------------------------------------
//

template <typename IntegralData>
FockMatrix<IntegralData>::FockMatrix(
      size_type n_alpha, size_type n_beta, const idata_type& integral_data,
      const linalgwrap::MultiVector<vector_type>& initial_guess, bool store_hf_terms)
      : m_n_alpha{n_alpha},
        m_n_beta{n_beta},
        m_idata_ptr{krims::make_subscription(integral_data, "FockMatrix")},
        m_fock_ptr(std::make_shared<stored_matrix_type>(integral_data.nbas(),
                                                        integral_data.nbas(), false)),
        m_store_hf_terms{store_hf_terms},
        m_terms_ptr{nullptr} {
  build_fock_matrix_from_coefficient(initial_guess);
}

template <typename IntegralData>
void FockMatrix<IntegralData>::calc_coulomb(const stored_matrix_type& density_bb,
                                            stored_matrix_type& coul_bb) const {
  assert_dbg(density_bb.n_rows() == density_bb.n_cols(), krims::ExcInternalError());
  assert_dbg(density_bb.n_rows() == coul_bb.n_rows(), krims::ExcInternalError());
  assert_dbg(density_bb.n_cols() == coul_bb.n_cols(), krims::ExcInternalError());

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

template <typename IntegralData>
void FockMatrix<IntegralData>::calc_exchange(const stored_matrix_type& density_bb,
                                             stored_matrix_type& exch_bb) const {
  assert_dbg(density_bb.n_rows() == density_bb.n_cols(), krims::ExcInternalError());

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

template <typename IntegralData>
void FockMatrix<IntegralData>::build_fock_matrix_from_density(stored_matrix_type pa_bb,
                                                              stored_matrix_type pb_bb) {
  // Assert densities have the correct size:
  assert_size(m_idata_ptr->nbas(), pa_bb.n_rows());
  assert_size(m_idata_ptr->nbas(), pa_bb.n_cols());
  assert_size(m_idata_ptr->nbas(), pb_bb.n_rows());
  assert_size(m_idata_ptr->nbas(), pb_bb.n_cols());

  if (!m_terms_ptr || !m_terms_ptr.unique()) {
    // Either the term storage is empty or we are not the only
    // one using it, so allocate new storage:
    m_terms_ptr = std::make_shared<terms_type>(m_idata_ptr->nbas(), false);
  }

  // Store densities:
  m_terms_ptr->pa_bb = std::move(pa_bb);
  m_terms_ptr->pb_bb = std::move(pb_bb);

  // Compute total density:
  stored_matrix_type pt_bb = m_terms_ptr->pa_bb + m_terms_ptr->pb_bb;

  // Build 2e terms:
  calc_coulomb(pt_bb, m_terms_ptr->j_bb);
  calc_exchange(m_terms_ptr->pa_bb, m_terms_ptr->ka_bb);
  calc_exchange(m_terms_ptr->pb_bb, m_terms_ptr->kb_bb);

  // Get shortcut references
  const stored_matrix_type& j_bb = m_terms_ptr->j_bb;
  const stored_matrix_type& ka_bb = m_terms_ptr->ka_bb;
  const stored_matrix_type& kb_bb = m_terms_ptr->kb_bb;

  // Calculate 1e energy:
  m_energies.energy_kinetic = trace(m_idata_ptr->t_bb() * pt_bb);
  m_energies.energy_elec_nuc_attr = trace(m_idata_ptr->v0_bb() * pt_bb);
  m_energies.energy_1e_terms =
        m_energies.energy_kinetic + m_energies.energy_elec_nuc_attr;

  // Calculate 2e energies and sum:
  m_energies.energy_coulomb = 0.5 * trace(j_bb * pt_bb);
  m_energies.energy_exchange = -0.5 * trace(ka_bb * m_terms_ptr->pa_bb);
  m_energies.energy_exchange -= 0.5 * trace(kb_bb * m_terms_ptr->pb_bb);
  m_energies.energy_2e_terms = m_energies.energy_coulomb + m_energies.energy_exchange;
  m_energies.energy_total = m_energies.energy_1e_terms + m_energies.energy_2e_terms;

  // If our fock pointer is not unique, i.e. there are other users,
  // make a new one first:
  if (!m_fock_ptr.unique()) {
    m_fock_ptr = std::make_shared<stored_matrix_type>(m_fock_ptr->n_rows(),
                                                      m_fock_ptr->n_cols(), false);
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

template <typename IntegralData>
void FockMatrix<IntegralData>::build_fock_matrix_from_coefficient(
      const linalgwrap::MultiVector<vector_type>& coefficients_bf) {
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
