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
#include "IntegralDataBase.hh"

namespace gscfmock {

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
  IntegralsSturmian14(scalar_type Z, scalar_type k_exp)
        : m_Z{Z},
          m_t_bb{(k_exp * k_exp) * m_t_bb_base},
          m_v0_bb{-k_exp * m_Z * m_v0_bb_base},
          m_i_bbbb{k_exp * m_i_bbbb_base} {}

  /** Return the number of basis functions */
  size_type nbas() const override { return m_nbas; }

  /** Return kinetic energy matrix */
  const matrix_type& t_bb() const override { return m_t_bb; }

  /** Return nuclear potential matrix */
  const matrix_type& v0_bb() const override { return m_v0_bb; }

  /** Return overlap matrix */
  const matrix_type& s_bb() const override { return m_s_bb; }

  /** Return two electron integrals as an nbas^2 x nbas^2 matrix
   *  of shell pairs.*/
  const matrix_type& i_bbbb() const override { return m_i_bbbb; }

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

}  // namespace gscfmock
