#pragma once
#include <linalgwrap/Exceptions.hh>
#include <string>

namespace gscf {

/** Class to represent the reason why an scf iteration failed
 *  Allows comparison with a ReasonId or another ScfFailReason
 *  in order to identify the reason from a caught exception */
class ScfFailReason {
public:
  enum class ReasonId {
    INNER_EIGENSOLVER_FAILED,
    MAXIMUM_ITERATIONS_REACHED,
    DIIS_STEP_FAILED,
  };

  /** Get the represented reason as a descriptive string */
  std::string as_string() const;

  /** Construct the class explicitly or implicitly from a
   * ReasonId */
  ScfFailReason(const ReasonId id) : m_id(id), m_further_info("") {}

  /** Construct the class explicitly or implicitly from a
   * ReasonId and provide some further information.*/
  ScfFailReason(const ReasonId id, const std::string further_info)
        : m_id(id), m_further_info(further_info) {}

  /** Compare the reason ids behind other and this */
  bool operator==(const ScfFailReason& other) const {
    return m_id == other.m_id;
  }

  /** Compare the reason ids behind other and this */
  bool operator!=(const ScfFailReason& other) const {
    return !this->operator==(other);
  }

private:
  const ReasonId m_id;
  const std::string m_further_info;
};

class ExcScfFailedToConverge : public ::linalgwrap::ExceptionBase {
public:
  ExcScfFailedToConverge(const ScfFailReason r) : m_fail_reason(r) {}
  virtual ~ExcScfFailedToConverge() noexcept {}

  virtual void print_extra(std::ostream& out) const noexcept;

private:
  const ScfFailReason m_fail_reason;
};

}  // namespace gscf
