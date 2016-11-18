#pragma once
#include <linalgwrap/Base/Solvers.hh>

namespace gscf {

/** Struct which contains the keys used for setting the
 *  ScfBase parameters */
struct ScfBaseKeys : public linalgwrap::IterativeSolverKeys {
  /** The number of eigenpairs to compute. Type: size_type */
  static const std::string n_eigenpairs;

  /** The submap containing the eigensolver parameters */
  static const std::string eigensolver_params;
};

}  // namespace gscf
