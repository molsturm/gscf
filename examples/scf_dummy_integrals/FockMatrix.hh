#pragma once
#include "integrals/IntegralDataBase.hh"
#include <linalgwrap/LazyMatrix_i.hh>
#include <linalgwrap/SubscriptionPointer.hh>
#include <linalgwrap/io.hh>
#include <linalgwrap/view.hh>

namespace scf_dummy {
// TODO more comments in the functions of this class

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

template <typename IntegralData>
class FockMatrix
      : public linalgwrap::LazyMatrix_i<typename IntegralData::matrix_type> {
    static_assert(std::is_base_of<IntegralDataBase, IntegralData>::value,
                  "IntegralData needs to be derived from IntegralDataBase");

  public:
    typedef IntegralData idata_type;
    typedef linalgwrap::LazyMatrix_i<typename idata_type::matrix_type>
          base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
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
     * TODO generalise.
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
     * 			       (exchange matrix, coulomb matrix, ...)
     */
    FockMatrix(size_type n_alpha, size_type n_beta,
               const idata_type& integral_data,
               const stored_matrix_type& initial_guess_bf,
               bool store_hf_terms = false);

    /** Return the number of rows of the matrix */
    size_type n_rows() const override;

    /** Return the number of columns of the matrix */
    size_type n_cols() const override;

    /** Return an element of the matrix */
    scalar_type operator()(size_type row, size_type col) const override;

    /** Apply matrix to another matrix */
    stored_matrix_type operator*(const stored_matrix_type& in) const override;

    /** Clone the matrix */
    lazy_matrix_expression_ptr_type clone() const override;

    /** Update the inner state:
     * Build the Fock matrix with the new coefficients
     *
     * It expects the new coefficients under the parameter key
     * "evec_coefficients"
     */
    void update(const linalgwrap::ParameterMap& map) override;

    /** Return the current set of energies, calculated at the most recent update
     */
    const energies_type& energies() const;

    /** Return the Hartree Fock terms
     *
     * asserts that they are actually stored in this class
     */
    const terms_type& hf_terms() const;

    /** Are the Hartree Fock terms stored? */
    bool are_hf_terms_stored() const;

  private:
    // Struct for the individual terms:

    /** Calculate the coulomb matrix from a given density */
    void calc_coulomb(const stored_matrix_type& density_bb,
                      stored_matrix_type& coul_bb) const;

    /** Calculate the exchange matrix from a given density */
    void calc_exchange(const stored_matrix_type& density_bb,
                       stored_matrix_type& exch_bb) const;

    /** Build the Fock matrix */
    void build_fock_matrix(const stored_matrix_type& coefficients_bf);

    /** Number of alpha electrons */
    size_type m_n_alpha;

    /** Number of beta electrons */
    size_type m_n_beta;

    //! The integral data
    linalgwrap::SubscriptionPointer<const idata_type> m_idata_ptr;

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
};

//
// --------------------------------------------------------------------
//

template <typename IntegralData>
FockMatrix<IntegralData>::FockMatrix(size_type n_alpha, size_type n_beta,
                                     const idata_type& integral_data,
                                     const stored_matrix_type& initial_guess,
                                     bool store_hf_terms)
      : m_n_alpha{n_alpha},
        m_n_beta{n_beta},
        m_idata_ptr{linalgwrap::make_subscription(integral_data, "FockMatrix")},
        m_fock_ptr(std::make_shared<stored_matrix_type>(
              integral_data.nbas(), integral_data.nbas(), false)),
        m_store_hf_terms{store_hf_terms},
        m_terms_ptr{nullptr} {
    build_fock_matrix(initial_guess);
}

template <typename IntegralData>
typename FockMatrix<IntegralData>::size_type FockMatrix<IntegralData>::n_rows()
      const {
    return m_fock_ptr->n_rows();
}

template <typename IntegralData>
typename FockMatrix<IntegralData>::size_type FockMatrix<IntegralData>::n_cols()
      const {
    return m_fock_ptr->n_cols();
}

template <typename IntegralData>
typename FockMatrix<IntegralData>::scalar_type FockMatrix<IntegralData>::
operator()(size_type row, size_type col) const {
    assert_greater(row, n_rows());
    assert_greater(col, n_cols());
    return (*m_fock_ptr)(row, col);
}

template <typename IntegralData>
typename FockMatrix<IntegralData>::stored_matrix_type FockMatrix<IntegralData>::
operator*(const stored_matrix_type& in) const {
    assert_size(m_fock_ptr->n_cols(), in.n_rows());
    return *m_fock_ptr * in;
}

template <typename IntegralData>
typename FockMatrix<IntegralData>::lazy_matrix_expression_ptr_type
FockMatrix<IntegralData>::clone() const {
    return lazy_matrix_expression_ptr_type(new FockMatrix<IntegralData>(*this));
}

template <typename IntegralData>
void FockMatrix<IntegralData>::update(const linalgwrap::ParameterMap& map) {
    // The coefficient key we look for:
    const std::string coeff_key = "evec_coefficients";

    if (map.exists(coeff_key)) {
        // We have new coefficients:
        const auto& coeff = map.at<stored_matrix_type>(coeff_key);
        build_fock_matrix(coeff);
    }
}

template <typename IntegralData>
inline const typename FockMatrix<IntegralData>::energies_type&
FockMatrix<IntegralData>::energies() const {
    return m_energies;
}

template <typename IntegralData>
inline const typename FockMatrix<IntegralData>::terms_type&
FockMatrix<IntegralData>::hf_terms() const {
    assert_dbg(m_terms_ptr, linalgwrap::ExcInvalidState("Terms not stored"));
    return *m_terms_ptr;
}

template <typename IntegralData>
inline bool FockMatrix<IntegralData>::are_hf_terms_stored() const {
    return m_store_hf_terms;
}

template <typename IntegralData>
void FockMatrix<IntegralData>::calc_coulomb(
      const stored_matrix_type& density_bb, stored_matrix_type& coul_bb) const {
    using namespace linalgwrap;

    assert_dbg(density_bb.n_rows() == density_bb.n_cols(), ExcInternalError());
    assert_dbg(density_bb.n_rows() == coul_bb.n_rows(), ExcInternalError());
    assert_dbg(density_bb.n_cols() == coul_bb.n_cols(), ExcInternalError());

    // Number of basis fuctions:
    size_t nbas = density_bb.n_rows();

    // J_{ab} = \sum_{cd} P_{cd} < ac | bd >
    // a and b are the same centre, so are c and d

    // Double loop over basis functions a and b:
    for (auto a : range(nbas)) {
        for (auto b : range(nbas)) {
            // Shell pair index for basis functions a and b:
            size_type ab_pair = a * nbas + b;

            // Sum accumulator variable:
            scalar_type sum{0};

            // Double loop over basis functions c and d:
            for (auto c : range(nbas)) {
                for (auto d : range(nbas)) {
                    // Shell pair index for basis functions c and d:
                    size_type cd_pair = c * nbas + d;

                    // Perform contraction:
                    sum += density_bb(c, d) *
                           m_idata_ptr->i_bbbb()(ab_pair, cd_pair);
                }  // d
            }      // c

            coul_bb(a, b) = sum;
        }  // b
    }      // a
}

template <typename IntegralData>
void FockMatrix<IntegralData>::calc_exchange(
      const stored_matrix_type& density_bb, stored_matrix_type& exch_bb) const {
    using namespace linalgwrap;

    assert_dbg(density_bb.n_rows() == density_bb.n_cols(), ExcInternalError());

    // Number of basis fuctions:
    size_t nbas = density_bb.n_rows();

    // Zero the exchange matrix:
    exch_bb.set_zero();

    // K_{ab} = \sum_{cd} P_{cd} < ab | cd >
    // a and c are the same centre, so are b and d

    // Double loop over shell pair basis functions a and c
    for (auto a : range(nbas)) {
        for (auto c : range(nbas)) {
            // Shell pair index for basis functions a and c:
            size_type ac_pair = a * nbas + c;

            // Double loop over shell pair basis functions d and b:
            for (auto d : range(nbas)) {
                for (auto b : range(nbas)) {
                    // Shell pair index for d and b:
                    size_type bd_pair = d * nbas + b;

                    // Perform contraction:
                    exch_bb(a, b) += density_bb(c, d) *
                                     m_idata_ptr->i_bbbb()(ac_pair, bd_pair);
                }  // b
            }      // d
        }          // c
    }              // a
}

template <typename IntegralData>
void FockMatrix<IntegralData>::build_fock_matrix(
      const stored_matrix_type& coefficients_bf) {
    using namespace linalgwrap;
    // Assert coefficients have the correct size:
    assert_size(m_idata_ptr->nbas(), coefficients_bf.n_rows());

    if (!m_terms_ptr || !m_terms_ptr.unique()) {
        // Either the term storage is empty or we are not the only
        // one using it, so allocate new storage:
        m_terms_ptr = std::make_shared<terms_type>(m_idata_ptr->nbas(), false);
    }

    // Coefficients for occupied orbitals only:
    auto ca_bo = view::columns(coefficients_bf, range(m_n_alpha));
    auto cb_bo = view::columns(coefficients_bf, range(m_n_beta));

    // Density for alpha and beta spin and the total density:
    m_terms_ptr->pa_bb =
          static_cast<stored_matrix_type>(ca_bo * view::transpose(ca_bo));
    m_terms_ptr->pb_bb =
          static_cast<stored_matrix_type>(cb_bo * view::transpose(cb_bo));
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
    // TODO this is a hell of a hack, due to the view::view call.
    // This does feel wrong. But otherwise we get a compile error!
    m_energies.energy_kinetic =
          (view::transpose(ca_bo) * view::view(m_idata_ptr->t_bb()) * ca_bo +
           view::transpose(cb_bo) * view::view(m_idata_ptr->t_bb()) * cb_bo)
                .trace();
    m_energies.energy_elec_nuc_attr =
          (view::transpose(ca_bo) * view::view(m_idata_ptr->v0_bb()) * ca_bo +
           view::transpose(cb_bo) * view::view(m_idata_ptr->v0_bb()) * cb_bo)
                .trace();
    m_energies.energy_1e_terms =
          m_energies.energy_kinetic + m_energies.energy_elec_nuc_attr;

    // Calculate 2e energies and sum:
    m_energies.energy_coulomb =
          0.5 * ((view::transpose(ca_bo) * view::view(j_bb) * ca_bo +
                  view::transpose(cb_bo) * view::view(j_bb) * cb_bo)
                       .trace());
    m_energies.energy_exchange =
          -0.5 * ((view::transpose(ca_bo) * view::view(ka_bb) * ca_bo +
                   view::transpose(cb_bo) * view::view(kb_bb) * cb_bo)
                        .trace());
    m_energies.energy_2e_terms =
          m_energies.energy_coulomb + m_energies.energy_exchange;
    m_energies.energy_total =
          m_energies.energy_1e_terms + m_energies.energy_2e_terms;

    // If our fock pointer is not unique, i.e. there are other users,
    // make a new one first:
    if (!m_fock_ptr.unique()) {
        m_fock_ptr = std::make_shared<stored_matrix_type>(
              m_fock_ptr->n_rows(), m_fock_ptr->n_cols(), false);
    }

    // Build two-electron terms and add everything to get fock matrix:
    // TODO adapt this method and class for open-shell calculations!
    *m_fock_ptr = m_idata_ptr->t_bb()  // kinetic energy term
                  +
                  m_idata_ptr->v0_bb()  // nuc-electron attract. potential term
                  + j_bb                // coloumb of full density
                  - ka_bb;              // exchange of same-spin density
    // TODO for beta fock matrix and uhf:
    //  m_fock_b = ... - build_exchange(pb_bb);

    // Deallocate storage if not needed any more:
    if (!m_store_hf_terms) {
        m_terms_ptr.reset();
    }
}

}  // namespace scf_dummy
