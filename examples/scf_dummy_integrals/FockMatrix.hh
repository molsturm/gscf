#pragma once
#include "integrals/IntegralDataBase.hh"
#include <linalgwrap/LazyMatrix_i.hh>
#include <linalgwrap/SubscriptionPointer.hh>
#include <linalgwrap/io.hh>
#include <linalgwrap/view.hh>  //XXX for debug

namespace scf_dummy {
// TODO more comments in the functions of this class

// XXX for debug
std::ofstream mathematicafile2("/tmp/debug_gscf_fock.m");
auto genoutfock = linalgwrap::io::make_writer<linalgwrap::io::Mathematica>(
      mathematicafile2, 1e-6);
// XXX for debug

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
     */
    FockMatrix(size_type n_alpha, size_type n_beta,
               const idata_type& integral_data,
               const stored_matrix_type& initial_guess_bf);

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

  private:
    /** Calculate the coulomb matrix from a given density */
    stored_matrix_type calc_coulomb(const stored_matrix_type& density_bb);

    /** Calculate the exchange matrix from a given density */
    stored_matrix_type calc_exchange(const stored_matrix_type& density_bb);

    /** Build the Fock matrix */
    void build_fock_matrix(const stored_matrix_type& coefficients_bf);

    /** Number of alpha electrons */
    size_type m_n_alpha;

    /** Number of beta electrons */
    size_type m_n_beta;

    //! The integral data
    linalgwrap::SubscriptionPointer<const idata_type> m_idata_ptr;

    //! The actual fock matrix for the current set of coefficients
    stored_matrix_type m_fock;
};

//
// --------------------------------------------------------------------
//

template <typename IntegralData>
FockMatrix<IntegralData>::FockMatrix(size_type n_alpha, size_type n_beta,
                                     const idata_type& integral_data,
                                     const stored_matrix_type& initial_guess)
      : m_n_alpha{n_alpha},
        m_n_beta{n_beta},
        m_idata_ptr{linalgwrap::make_subscription(integral_data, "FockMatrix")},
        m_fock(integral_data.nbas(), integral_data.nbas(), false) {
    build_fock_matrix(initial_guess);
}

template <typename IntegralData>
typename FockMatrix<IntegralData>::size_type FockMatrix<IntegralData>::n_rows()
      const {
    return m_fock.n_rows();
}

template <typename IntegralData>
typename FockMatrix<IntegralData>::size_type FockMatrix<IntegralData>::n_cols()
      const {
    return m_fock.n_cols();
}

template <typename IntegralData>
typename FockMatrix<IntegralData>::scalar_type FockMatrix<IntegralData>::
operator()(size_type row, size_type col) const {
    assert_greater(row, n_rows());
    assert_greater(col, n_cols());
    return m_fock(row, col);
}

template <typename IntegralData>
typename FockMatrix<IntegralData>::stored_matrix_type FockMatrix<IntegralData>::
operator*(const stored_matrix_type& in) const {
    assert_size(m_fock.n_cols(), in.n_rows());
    return m_fock * in;
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
typename FockMatrix<IntegralData>::stored_matrix_type
FockMatrix<IntegralData>::calc_coulomb(const stored_matrix_type& density_bb) {
    using namespace linalgwrap;

    assert_dbg(density_bb.n_rows() == density_bb.n_cols(), ExcInternalError());

    // Number of basis fuctions:
    size_t nbas = density_bb.n_rows();

    // Output coulomb matrix:
    stored_matrix_type j_bb(nbas, nbas, false);

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

            j_bb(a, b) = sum;
        }  // b
    }      // a

    return j_bb;
}

template <typename IntegralData>
typename FockMatrix<IntegralData>::stored_matrix_type
FockMatrix<IntegralData>::calc_exchange(const stored_matrix_type& density_bb) {
    using namespace linalgwrap;

    // TODO something not quite correct here ...

    assert_dbg(density_bb.n_rows() == density_bb.n_cols(), ExcInternalError());

    // Number of basis fuctions:
    size_t nbas = density_bb.n_rows();

    // Output exchange matrix:
    stored_matrix_type k_bb(nbas, nbas, true);

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
                    k_bb(a, b) += density_bb(c, d) *
                                  m_idata_ptr->i_bbbb()(ac_pair, bd_pair);
                }  // b
            }      // d
        }          // c
    }              // a

    return k_bb;
}

template <typename IntegralData>
void FockMatrix<IntegralData>::build_fock_matrix(
      const stored_matrix_type& coefficients_bf) {
    using namespace linalgwrap;

    // Assert coefficients have the correct size:
    assert_size(m_idata_ptr->nbas(), coefficients_bf.n_rows());

    // Coefficients for occupied orbitals only:
    auto occa = view::columns(coefficients_bf, range(m_n_alpha));
    auto occb = view::columns(coefficients_bf, range(m_n_beta));

    // Resulting density for alpha and beta spin:
    stored_matrix_type pa_bb =
          static_cast<stored_matrix_type>(occa * view::transpose(occa));
    stored_matrix_type pb_bb =
          static_cast<stored_matrix_type>(occb * view::transpose(occb));

    // Total density
    stored_matrix_type pt_bb = pa_bb + pb_bb;

    // Build two-electron terms and add everything to get fock matrix:

    // TODO adapt this method and class for open-shell calculations!
    m_fock = m_idata_ptr->t_bb()      // kinetic energy term
             + m_idata_ptr->v0_bb()   // nuc-electron attract. potential term
             + calc_coulomb(pt_bb)    // coloumb of full density
             - calc_exchange(pa_bb);  // exchange of same-spin density
    // TODO for beta fock matrix and uhf:
    //  m_fock_b = ... - build_exchange(pb_bb);

    // XXX for debug
    static size_t count = 0;
    genoutfock.write(coefficients_bf, "cbf" + std::to_string(++count));
    genoutfock.write(pt_bb, "pbb" + std::to_string(count));
    genoutfock.write(calc_coulomb(pt_bb), "jbb" + std::to_string(count));
    genoutfock.write(calc_exchange(pa_bb), "kbb" + std::to_string(count));
    genoutfock.write(m_fock, "fbb" + std::to_string(count));
    // XXX for debug
}

}  // namespace scf_dummy
