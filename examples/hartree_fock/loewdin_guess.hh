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
#include <linalgwrap/MultiVector.hh>

namespace dummy_scf {
using gscfmock::matrix_type;
using gscfmock::scalar_type;
using gscfmock::vector_type;

linalgwrap::MultiVector<vector_type> loewdin_guess(const matrix_type& overlap_bb);

}  // namespace dummy_scf