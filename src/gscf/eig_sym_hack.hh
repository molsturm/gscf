#pragma once

#include <linalgwrap/ArmadilloMatrix.hh>
#include <linalgwrap/view.hh>

namespace gscf {

namespace detail {

// XXX: Hack an inverse square root function until its in linalgwrap
template <typename Matrix>
Matrix invsqrtmat_hack(const Matrix& m) {
  using namespace linalgwrap;

  // Assert some types:
  static_assert(std::is_same<Matrix, ArmadilloMatrix<double>>::value,
                "This hack only works for armadillo matrices");

  // Assert that m is symmetric
  assert_dbg(m.is_symmetric(), ExcMatrixNotSymmetric());

  // Unpack matrix:
  arma::mat m_arma = m.data().t();

  // Diagonalise:
  arma::vec eval_arma;
  arma::mat evec_arma;
  bool success = arma::eig_sym(eval_arma, evec_arma, m_arma);
  assert_abort(success, ExcInternalError());

#ifdef DEBUG
  // Assert that the input matrix is positive definite:
  for (auto diag : eval_arma) {
    assert_greater_equal(0, diag);
  }
#endif

  // take 1/sqrt( . ) for each eigenvalue:
  std::transform(std::begin(eval_arma), std::end(eval_arma),
                 std::begin(eval_arma),
                 [](double elem) { return 1. / std::sqrt(elem); });

  // Construct new matrix:
  arma::mat res_arma = evec_arma * arma::diagmat(eval_arma) * evec_arma.t();

  // Enwrap (assuming that res_arma is symmetric, which it is by construction)
  return ArmadilloMatrix<double>(std::move(res_arma));
}

// XXX: Hack an eigensolver until its in linalgwrap
template <typename ProblemMatrix>
bool eig_sym_hack(
      const ProblemMatrix& f,
      const typename ProblemMatrix::stored_matrix_type& metricmat_bb,
      typename ProblemMatrix::stored_matrix_type& evec,
      linalgwrap::VectorOf<typename ProblemMatrix::stored_matrix_type>& eval) {
  using namespace linalgwrap;

  assert_dbg(f.is_symmetric(), ExcMatrixNotSymmetric());
  assert_dbg(metricmat_bb.is_symmetric(), ExcMatrixNotSymmetric());

  // Assert some types:
  static_assert(std::is_same<typename ProblemMatrix::stored_matrix_type,
                             ArmadilloMatrix<double>>::value,
                "This hack only works for armadillo matrices");

  // Take inverse square root of metric:
  ArmadilloMatrix<double> invsqrt_metric = invsqrtmat_hack(metricmat_bb);

  // Assert that this is symmetric:
  assert_dbg(invsqrt_metric.is_symmetric(), ExcMatrixNotSymmetric());

  // Transform the fock matrix into the orthogonal basis:
  // Assume that invsqrt_metric is hermetian (as it should)
  auto fock_transformed = invsqrt_metric * (f * invsqrt_metric);

  // Get the transformed fock matrix in memory:
  ArmadilloMatrix<double> t_fock_mem =
        static_cast<ArmadilloMatrix<double>>(fock_transformed);

  // Unpack matrices:
  // Note that the .t() is skipped since fock is symmetric
  arma::mat fock_arma = t_fock_mem.data();

  // Diagonalise:
  arma::vec eval_arma;  // This is a column matrix
  arma::mat evec_arma;
  bool result = arma::eig_sym(eval_arma, evec_arma, fock_arma);
  if (!result) {
    // arma encountered an error
    return false;
  }

  // Properly repack and transform back to original basis:
  evec = invsqrt_metric * ArmadilloMatrix<double>(evec_arma.t());
  eval = VectorOf<ArmadilloMatrix<double>>(std::begin(eval_arma),
                                           std::end(eval_arma));

  return true;
}

}  // namespace detail

}  // namespace gscf
