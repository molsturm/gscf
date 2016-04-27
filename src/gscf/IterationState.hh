#pragma once
#include <cstddef>
#include <string>
namespace gscf {

class IterationState {
  public:
    typedef size_t count_type;

    /** \brief Default constructor
     *
     * The failbit is unset and the iteration count is 0.
     * */
    IterationState();

    /** Fail the iteraton and specify a reason why */
    void fail(std::string reason);

    /** Clear the failed flag and clear the fail reason. */
    void clear_failed();

    /** Increase the internal iteration count by one
     *  and return the new iteration count. */
    count_type increase_iteration_count();

    /** Is the failed bit set? */
    bool is_failed() const;

    /** Get the fail reason
     *
     * \note Returns an empty string in case that the iteration has not failed.
     * An empty string does, however, not imply that the iteration has not
     * failed.
     */
    const std::string& fail_reason() const;

    /** Return the current iteration count */
    count_type n_iter_count() const;

  private:
    count_type m_count;

    //! Has the iteration failed?
    bool m_failed;

    //! Why has it failed?
    std::string m_fail_reason;
};

//
// ---------------------------------------------------------
//

inline IterationState::IterationState()
      : m_count{0}, m_failed{false}, m_fail_reason{} {}

inline void IterationState::fail(std::string reason) {
    m_failed = true;
    m_fail_reason = reason;
}

inline void IterationState::clear_failed() {
    m_failed = false;
    m_fail_reason = "";
}

inline const std::string& IterationState::fail_reason() const {
    return m_fail_reason;
}

inline typename IterationState::count_type
IterationState::increase_iteration_count() {
    return ++m_count;
}

inline bool IterationState::is_failed() const { return m_failed; }

inline typename IterationState::count_type IterationState::n_iter_count()
      const {
    return m_count;
}

}  // namespace gscf
