#pragma once
#include "ScfBaseKeys.hh"

namespace gscf {

/** Struct which contains the keys used for setting the
 *  ScfBase parameters */
struct PulayDiisScfKeys : public ScfBaseKeys {
  /** The number of previous SCF steps to consider for the DIIS.
   *  Type: size_type */
  static const std::string n_prev_steps;

  /** The DIIS error which is allowed to reach convergence
   *  Type: real_type */
  static const std::string max_error_norm;
};

}  // namespace gscf
