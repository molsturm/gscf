#pragma once
#include <linalgwrap/ArmadilloMatrix.hh>
#include <linalgwrap/VectorOf.hh>

namespace gscf {
namespace detail {

// XXX: Hack a linear solver until its in linalgwrap
template <typename StoredMatrix>
bool linsolve_hack(const StoredMatrix& m,
                   const linalgwrap::VectorOf<StoredMatrix>& rhs,
                   linalgwrap::VectorOf<StoredMatrix>& x) {
  using namespace linalgwrap;

  assert_dbg(m.is_symmetric(), ExcMatrixNotSymmetric());
  assert_size(m.n_rows(), rhs.size());
  assert_size(m.n_cols(), x.size());

  // Assert some types:
  static_assert(std::is_same<StoredMatrix, ArmadilloMatrix<double>>::value,
                "This hack only works for armadillo matrices");

  // Unpack matrices:
  arma::mat m_arma = m.data();  //.t() is skipped since m is symmetric
  arma::vec rhs_arma(rhs.data().t());

  // Solve the linear system:
  arma::vec x_arma;
  bool result = arma::solve(x_arma, m_arma, rhs_arma);
  if (!result) {
    // arma encountered an error
    return false;
  }

  // Properly repack:
  x = VectorOf<ArmadilloMatrix<double>>(std::begin(x_arma), std::end(x_arma));
  return true;
}

}  // namespace detail
}  // namespace gscf
