#pragma once
#include "IntegralDataBase.hh"

namespace scf_dummy {

/** \brief Dummy integral class with 14 sturmians
 *
 * Corresponds to an atomic Sturmian Basis set with
 * a maximum principle quantum number of 3.
 */
class IntegralsSturmian14 : public IntegralDataBase {
  public:
    typedef IntegralDataBase base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::matrix_type matrix_type;
    typedef typename base_type::size_type size_type;

    /** \name Construct the data object
     *
     * \param Z      nuclear charge
     * \param k_exp  Sturmian exponent factor
     * */
    IntegralsSturmian14(scalar_type Z, scalar_type k_exp);

    /** Return the number of basis functions */
    size_type nbas() const override;

    /** Return kinetic energy matrix */
    const matrix_type& t_bb() const override;

    /** Return nuclear potential matrix */
    const matrix_type& v0_bb() const override;

    /** Return overlap matrix */
    const matrix_type& s_bb() const override;

    /** Return two electron integrals as an nbas^2 x nbas^2 matrix
     *  of shell pairs.*/
    const matrix_type& i_bbbb() const override;

  private:
    static constexpr size_type m_nbas = 14;
    const scalar_type m_Z;       //< number of nuclei
    const matrix_type m_t_bb;    //< Actual t_bb matrix to be returned.
    const matrix_type m_v0_bb;   //< Actual v0_bb matrix to be returned.
    const matrix_type m_i_bbbb;  //< Actual i_bbbb matrix to be returned.

    //! Kinetic matrix without the k or Z factors.
    static const matrix_type m_t_bb_base;

    //! Overlap matrix
    static const matrix_type m_s_bb;

    //! Electron-core interaction matrix without k or Z factors.
    static const matrix_type m_v0_bb_base;

    //! Two electron integrals without k or Z factors.
    static const matrix_type m_i_bbbb_base;
};

}  // namespace scf_dummy
