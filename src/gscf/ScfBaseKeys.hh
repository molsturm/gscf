#pragma once
#include <linalgwrap/Base/Solvers.hh>

namespace gscf {

/** Struct which contains the keys used for setting the
 *  ScfBase parameters */
struct ScfBaseKeys : public linalgwrap::IterativeWrapperKeys {
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
