#pragma once
#include "IterationState.hh"
#include <linalgwrap/Constants.hh>

namespace gscf {

template <typename SizeType>
struct IterationConstants {
  /** Number representing the notion, that all states, eigenvalues or whatever
   *  are to be seeked by the algorithm */
  static constexpr SizeType all = linalgwrap::Constants<SizeType>::invalid;
};
}
