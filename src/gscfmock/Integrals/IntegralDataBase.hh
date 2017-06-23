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
#include <gscfmock/types.hh>

namespace gscfmock {

/** Class for dummy integrals around an atom */
class IntegralDataBase : public krims::Subscribable {
 public:
  IntegralDataBase()                        = default;
  virtual ~IntegralDataBase()               = default;
  IntegralDataBase(IntegralDataBase&&)      = default;
  IntegralDataBase(const IntegralDataBase&) = default;
  IntegralDataBase& operator=(IntegralDataBase&&) = default;
  IntegralDataBase& operator=(const IntegralDataBase&) = default;

  virtual size_t nbas() const = 0;

  /** \name Data access */
  ///@{
  /** Return kinetic energy matrix */
  virtual const matrix_type& t_bb() const = 0;

  /** Return nuclear potential matrix */
  virtual const matrix_type& v0_bb() const = 0;

  /** Return overlap matrix */
  virtual const matrix_type& s_bb() const = 0;

  /** Return two electron integrals as an nbas^2 x nbas^2 matrix
   *  of shell pairs.*/
  virtual const matrix_type& i_bbbb() const = 0;
};
}  // namespace gscfmock
