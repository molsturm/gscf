#include "ExcScfFailedToConverge.hh"

namespace gscf {

std::string ScfFailReason::as_string() const {
  std::string base;
  switch (m_id) {
    case ReasonId::INNER_EIGENSOLVER_FAILED:
      base = "Inner eigenproblem failed to converge";
      break;
    case ReasonId::MAXIMUM_ITERATIONS_REACHED:
      base = "Maximum number of iterations reached.";
      break;
    case ReasonId::DIIS_STEP_FAILED:
      base = "The DIIS accellerator step has failed.";
      break;
  }

  return base + m_further_info;
}

void ExcScfFailedToConverge::print_extra(std::ostream& out) const noexcept {
  out << "The SCF procedure failed to converge: " << m_fail_reason.as_string()
      << std::endl;
}

}  // namespace gscf
