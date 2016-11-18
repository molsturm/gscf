#pragma once
#include <linalgwrap/TypeUtils.hh>
#include <linalgwrap/eigensystem.hh>

namespace scf_dummy {
using namespace linalgwrap;

template <typename Matrix>
MultiVector<typename Matrix::vector_type> loewdin_guess(
      const Matrix& overlap_bb) {
  typedef typename Matrix::size_type size_type;
  typedef typename Matrix::vector_type vector_type;

  // apply Löwdin normalisation to the basis functions
  //   - Diagonalise the overlap
  //   - Take 1/\sqrt{evals} at the diagonal
  //   - results in orthonormalised basis functions

  try {
    auto sol = eigensystem_hermitian(overlap_bb);

    // Eigenvectors and eigenvalues.
    auto& evectors = sol.evectors();
    const auto& evalues = sol.evalues();

    assert_dbg(evectors.n_vectors(), overlap_bb.n_cols());
    assert_dbg(evectors.n_elem(), overlap_bb.n_rows());

    for (size_type i = 0; i < evectors.n_vectors(); ++i) {
      evectors[i] *= 1. / sqrt(evalues[i]);
    }

    return std::move(sol.evectors());
  } catch (const linalgwrap::SolverException& e) {
    std::cerr << "Obtaining Löwdin guess failed, using zero guess" << std::endl;
    return MultiVector<vector_type>(overlap_bb.n_rows(), overlap_bb.n_cols());
  }
}

}  // namespace scf_dummy
