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
#include <linalgwrap/BlockDiagonalMatrix.hh>
#include <linalgwrap/LazyMatrixSum.hh>

namespace gscf {

template <typename Operator,
          bool BlockDiagonal = linalgwrap::IsBlockDiagonalMatrix<Operator>::value>
struct OperatorLinearCombination
      : public linalgwrap::LazyMatrixSum<typename Operator::stored_matrix_type> {
  typedef Operator operator_type;
  typedef typename operator_type::scalar_type scalar_type;
  typedef linalgwrap::LazyMatrixSum<typename Operator::stored_matrix_type> base_type;

  /** \note An operator is required at initialisation since the BlockDiagonalMatrix
   *        need to know the sizes of the blocks
   */
  explicit OperatorLinearCombination(const operator_type& op, scalar_type coeff = 1)
        : base_type{op, coeff} {}

  /**  Add a term to the linear combination
   *
   * It is assumed that the terms match in size
   */
  void push_term(const operator_type& op, scalar_type coeff = 1) {
    (*this) += coeff * op;
  }
};

template <typename Operator>
struct OperatorLinearCombination<Operator, true>
      : public linalgwrap::BlockDiagonalMatrix<
              linalgwrap::LazyMatrixSum<typename Operator::stored_matrix_type>,
              Operator::n_blocks> {
  typedef Operator operator_type;
  typedef typename operator_type::scalar_type scalar_type;
  typedef linalgwrap::BlockDiagonalMatrix<
        linalgwrap::LazyMatrixSum<typename Operator::stored_matrix_type>,
        Operator::n_blocks>
        base_type;
  typedef linalgwrap::LazyMatrixSum<typename Operator::stored_matrix_type> block_type;

  /** \note An operator is required at initialisation since the BlockDiagonalMatrix
   *        need to know the sizes of the blocks
   */
  explicit OperatorLinearCombination(const operator_type& op, scalar_type coeff = 1)
        : base_type{{{block_type{op.diag_blocks()[0], coeff},
                      block_type{op.diag_blocks()[1], coeff}}}} {
    static_assert(Operator::n_blocks == 2,
                  "This version of OperatorLinearCombination assumes that there are "
                  "exactly two block.");
  }

  /**  Add a term to the linear combination
   *
   * It is assumed that the terms match in size
   */
  void push_term(const operator_type& op, scalar_type coeff = 1) {
    auto itop   = std::begin(op.diag_blocks());
    auto itthis = std::begin(this->diag_blocks());

    for (; itop != std::end(op.diag_blocks()); ++itop, ++itthis) {
      itthis->push_term(*itop, coeff);
    }
  }
};

}  // namespace gscf
