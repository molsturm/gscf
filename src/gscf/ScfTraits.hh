#pragma once
#include <linalgwrap/VectorOf.hh>

namespace gscf {

/** Default traits for the scf classes */
template <typename ProblemMatrix>
struct ScfTraits {
  typedef ProblemMatrix probmat_type;
  typedef typename probmat_type::scalar_type scalar_type;
  typedef typename probmat_type::size_type size_type;
  typedef typename probmat_type::stored_matrix_type matrix_type;
  typedef linalgwrap::VectorOf<matrix_type> vector_type;
};

}  // namespace gscf
