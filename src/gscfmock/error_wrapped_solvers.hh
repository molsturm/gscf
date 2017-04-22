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
#include "FockMatrix.hh"
#include "pulay_error.hh"
#include <gscf/PlainScf.hh>
#include <gscf/PulayDiisScf.hh>
#include <gscf/TruncatedOptDampScf.hh>

namespace gscfmock {
namespace error_wrapped_solvers {

template <typename Solver>
class ErrorWrapped : public Solver {
 public:
  using Solver::Solver;

  matrix_type calculate_error(const typename Solver::state_type& s) const override {
    return pulay_error(s.problem_matrix(), s.eigensolution().evectors(),
                       s.overlap_matrix());
  }
};

using PlainScf =
      ErrorWrapped<gscf::PlainScf<gscf::PlainScfState<FockMatrix, matrix_type>>>;
using PulayDiisScf =
      ErrorWrapped<gscf::PulayDiisScf<gscf::PulayDiisScfState<FockMatrix, matrix_type>>>;
using TruncatedOptDampScf = ErrorWrapped<
      gscf::TruncatedOptDampScf<gscf::TruncatedOptDampScfState<FockMatrix, matrix_type>>>;
}  // namespace error_wrapped_solvers
}  // namespace gscfmock
