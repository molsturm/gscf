#pragma once
#include "IterationControl.hh"

namespace gscf {

/** Class to provide the basic set of paramters the ScfBase class needs.
 *
 * \see gscf::ScfBase
 */
template <typename ScfTraits>
struct ScfControlBase : public IterationControl {
    typedef ScfTraits scf_traits_type;
    typedef typename scf_traits_type::probmat_type probmat_type;
    typedef typename scf_traits_type::scalar_type scalar_type;
    typedef typename scf_traits_type::size_type size_type;
    typedef typename scf_traits_type::matrix_type matrix_type;
    typedef typename scf_traits_type::vector_type vector_type;

    /** The number of eigenpairs to calculate in the SCF procedure */
    size_type n_eigenpairs;

    /** Construct an ScfControlBase object and set some defaults
     *
     * \param n_eigenpairs_    Number of eigenpairs to compute
     *                         in the SCF procedure
     * \param max_iter_        The maximum number of iterations to perform.
     * */
    ScfControlBase(size_type n_eigenpairs_, size_type max_iter_)
          : IterationControl(max_iter_), n_eigenpairs(n_eigenpairs_){};
};

}  // gscf
