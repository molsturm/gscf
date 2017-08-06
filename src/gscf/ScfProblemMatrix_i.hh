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
#include <string>

namespace gscf {

/** This class defines the extra interface, which all SCF solvers
 *  expect a problem matrix to satisfy.
 *
 *  Apart from the interface given here the problem matrix should
 *  be a lazyten::LazyMatrixExpression, i.e. its internal state
 *  should be updateable using the function update().
 */
class ScfProblemMatrix_i {
 public:
  /** Return the key which should be used by the SCF solver to update the state of the
   * Problem matrix with the new eigenvector coefficients */
  virtual const std::string& scf_update_key() const = 0;

  ScfProblemMatrix_i()                          = default;
  virtual ~ScfProblemMatrix_i()                 = default;
  ScfProblemMatrix_i(const ScfProblemMatrix_i&) = default;
  ScfProblemMatrix_i(ScfProblemMatrix_i&&)      = default;
  ScfProblemMatrix_i& operator=(ScfProblemMatrix_i&&) = default;
  ScfProblemMatrix_i& operator=(const ScfProblemMatrix_i&) = default;
};

}  // namespace gscf
