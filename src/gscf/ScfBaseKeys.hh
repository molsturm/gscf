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
#include <lazyten/Base/Solvers.hh>

namespace gscf {

/** Struct which contains the keys used for setting the
 *  ScfBase parameters */
struct ScfBaseKeys : public lazyten::IterativeWrapperKeys {
  /** The number of eigenpairs to compute. Type: size_type */
  static const std::string n_eigenpairs;

  /** The submap containing the eigensolver parameters */
  static const std::string eigensolver_params;

  /** The norm of the SCF error matrix/vector which is allowed
   *  in order to reach convergence.
   *  Type: real_type */
  static const std::string max_error_norm;
};

}  // namespace gscf
