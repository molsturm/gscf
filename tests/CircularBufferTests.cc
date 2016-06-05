#include <catch.hpp>
#include <gscf/CircularBuffer.hh>
#include <list>
#include <rapidcheck.h>
#include <rapidcheck/state.h>

// have an extra verbose output for rapidcheck function tests:
//#define HAVE_CIRCBUFF_RC_CLASSIFY

namespace gscf {
namespace tests {
using namespace rc;

namespace circular_buffer_tests {
/* Model for testing CircularBuffer */
template <typename T>
struct CircularBufferModel {
  /** Maximum size of the circular buffer */
  size_t max_size;

  /** List containing the actual data
   * in the actual expected order.
   * Contains size elements */
  std::list<T> data;

  void assert_equivalent_to(const CircularBuffer<T>& b) const {
    RC_ASSERT(b.max_size() == max_size);
    RC_ASSERT(b.size() == data.size());
    RC_ASSERT(b.empty() == data.empty());

    if (b.size() > 0) {
      RC_ASSERT(b.front() == data.front());
      RC_ASSERT(b.back() == data.back());

      auto itdata = std::begin(data);
      for (auto it = b.begin(); it != b.end(); ++it, ++itdata) {
        RC_ASSERT(*itdata == *it);
      }
    }  // if
  }
};

//
// Operations and commands
//

/** Push a random object to the front of the buffer */
template <typename T>
struct PushFront
      : rc::state::Command<CircularBufferModel<T>, CircularBuffer<T>> {
  typedef CircularBufferModel<T> model_type;
  typedef CircularBuffer<T> sut_type;

  T t;  //< New object to push_front;
  PushFront() : t{*gen::arbitrary<T>()} {};

  void apply(model_type& model) const override {
    RC_PRE(model.max_size > 0u);
    model.data.push_front(t);

    if (model.data.size() > model.max_size) {
      model.data.resize(model.max_size);
    }
  }

  void run(const model_type& model, sut_type& sut) const override {
#ifdef HAVE_CIRCBUFF_RC_CLASSIFY
    RC_CLASSIFY(true, "PushFront");
#endif

    // Perform on sut
    sut.push_front(t);

    // Check with model:
    auto nmodel = this->nextState(model);
    nmodel.assert_equivalent_to(sut);
  }

  void show(std::ostream& os) const override {
    os << "PushFront (" << t << ")";
  }
};  // PushFront

/** Push a random object to the back of the buffer */
template <typename T>
struct PushBack
      : rc::state::Command<CircularBufferModel<T>, CircularBuffer<T>> {
  typedef CircularBufferModel<T> model_type;
  typedef CircularBuffer<T> sut_type;

  T t;  //< New object to push_back;
  PushBack() : t{*gen::arbitrary<T>()} {};

  void apply(model_type& model) const override {
    RC_PRE(model.max_size > 0u);
    model.data.push_back(t);

    while (model.data.size() > model.max_size) {
      model.data.pop_front();
    }
  }

  void run(const model_type& model, sut_type& sut) const override {
#ifdef HAVE_CIRCBUFF_RC_CLASSIFY
    RC_CLASSIFY(true, "PushBack");
#endif

    // Perform on sut
    sut.push_back(t);

    // Check with model:
    auto nmodel = this->nextState(model);
    nmodel.assert_equivalent_to(sut);
  }

  void show(std::ostream& os) const override { os << "PushBack (" << t << ")"; }
};  // PushBack

/** Clear the buffer */
template <typename T>
struct Clear : rc::state::Command<CircularBufferModel<T>, CircularBuffer<T>> {
  typedef CircularBufferModel<T> model_type;
  typedef CircularBuffer<T> sut_type;

  void apply(model_type& model) const override {
    RC_PRE(model.max_size > 0u);
    model.data.clear();
  }

  void run(const model_type& model, sut_type& sut) const override {
#ifdef HAVE_CIRCBUFF_RC_CLASSIFY
    RC_CLASSIFY(true, "Clear");
#endif

    // Perform on sut
    sut.clear();

    // Check with model:
    this->nextState(model).assert_equivalent_to(sut);
  }

  void show(std::ostream& os) const override { os << "Clear"; }
};  // Clear

/** Change max_size */
template <typename T>
struct ChangeMaxSize
      : rc::state::Command<CircularBufferModel<T>, CircularBuffer<T>> {
  typedef CircularBufferModel<T> model_type;
  typedef CircularBuffer<T> sut_type;

  size_t max;  //< New max_size
  ChangeMaxSize() : max{*gen::inRange<size_t>(0, 11)} {}

  void apply(model_type& model) const override {
    model.max_size = max;
    if (model.data.size() > max) {
      model.data.resize(max);
    }
  }

  void run(const model_type& model, sut_type& sut) const override {
#ifdef HAVE_CIRCBUFF_RC_CLASSIFY
    RC_CLASSIFY(true, "ChangeMaxSize" + std::to_string(max));
#endif

    // Perform on sut
    sut.max_size(max);

    // Check with model:
    this->nextState(model).assert_equivalent_to(sut);
  }

  void show(std::ostream& os) const override {
    os << "ChangeMaxSize (" << max << ")";
  }
};  // ChangeMaxSize

/** execute a random test of a list of the above commands
 *
 * \tparam T typename of the circular buffer
 * \tparam Commands The command list
 */
template <typename T, typename... Commands>
void exectute_random_test() {
  size_t max_size = *rc::gen::inRange<size_t>(2, 10).as("max_size");
  std::list<T> data{};

  // Setup the model and initial system state
  typedef CircularBufferModel<T> model_type;
  typedef CircularBuffer<T> sut_type;
  model_type model{};
  model.max_size = max_size;
  model.data = data;
  sut_type sut{max_size};

  // Generate command sequence and execute random state test
  auto genCommands = rc::state::gen::execOneOfWithArgs<Commands...>;
  state::check(model, sut, genCommands());
};

}  // namespace circular_buffer_tests

//
// ---------------------------------------------------------------
//

template <typename Iterator>
void showValue(const gscf::CircularIterator<Iterator>& it, std::ostream& os) {
  os << "Position element: " << *it.position() << std::endl
     << "Begin element: " << *it.iteration_range().first << std::endl;
}

//
// ---------------------------------------------------------------
//

TEST_CASE("Circular Iterator and Circular Buffer", "[circular_buffer]") {
  using namespace circular_buffer_tests;

  // Make sure that the program does not get aborted
  linalgwrap::AssertDbgEffect::set(linalgwrap::ExceptionEffect::THROW);

  // The type to do the tests with
  typedef int testtype;

  SECTION("Incrementing and decrementing an empty range") {
    auto test = [](std::vector<testtype> v) {
      size_t pos = 0;
      if (v.size() > 0) {
        pos = *gen::inRange<size_t>(0, v.size())
                     .as("position of the empty range.");
      }

      // Construct the empty range iterator:
      typedef CircularIterator<decltype(v.begin())> circit_t;
      circit_t it{std::begin(v) + pos, std::begin(v) + pos,
                  std::begin(v) + pos};

#ifdef DEBUG
      RC_ASSERT_THROWS_AS(++it, linalgwrap::ExcInvalidState);
      RC_ASSERT_THROWS_AS(--it, linalgwrap::ExcInvalidState);
      RC_ASSERT_THROWS_AS(it++, linalgwrap::ExcInvalidState);
      RC_ASSERT_THROWS_AS(it--, linalgwrap::ExcInvalidState);
      RC_ASSERT_THROWS_AS(*it, linalgwrap::ExcInvalidState);
#endif
    };

    REQUIRE(rc::check(
          "CircularIterator: Incrementing and decrementing an empty range",
          test));
  }

  //
  // ---------------------------------------------------------------
  //

  SECTION("Check that pre- and post-decrement/increment agree") {
    auto test = [](std::vector<testtype> v) {
      RC_PRE(v.size() > 0u);
      size_t pos = *gen::inRange<size_t>(0, v.size())
                          .as("position to start the circle");

      auto it = circular_begin(v, pos);
      auto itp1 = ++it;
      auto itp2 = it++;
      RC_ASSERT(itp1 == itp2);

      auto itm1 = --it;
      auto itm2 = it--;
      RC_ASSERT(itm1 == itm2);
    };

    REQUIRE(rc::check(
          "CircularIterator: Agreement of pre-/post- decrement/increment",
          test));
  }

  //
  // ---------------------------------------------------------------
  //

  SECTION("Check that ++(--it) and --(++it)) are identical.") {
    auto test = [](std::vector<testtype> v) {
      RC_PRE(v.size() > 0u);
      size_t pos = *gen::inRange<size_t>(0, v.size())
                          .as("position to start the circle");

      auto it = circular_begin(v, pos);
      auto res1 = --(++it);
      RC_ASSERT(res1 == it);

      auto res2 = ++(--it);
      RC_ASSERT(res2 == it);

      RC_ASSERT(*res1 == *res2);
      RC_ASSERT(res1 == res2);
    };

    REQUIRE(rc::check(
          "CircularIterator: Check that ++(--it) and --(++it)) are identical.",
          test));
  }

  //
  // ---------------------------------------------------------------
  //

  SECTION("Infinite circular iteration around arbitrary vector.") {
    auto test = [](std::vector<testtype> v) {
      RC_PRE(v.size() > 0u);
      RC_PRE(v.size() < 11u);

      size_t iterations = *gen::inRange<size_t>(1, 10 * v.size())
                                 .as("Number of iterations to carry out");
      size_t startpos = *gen::inRange<size_t>(0, v.size())
                               .as("position to start the circle");

      // Test forward
      auto it = circular_begin(v, startpos);
      for (size_t i = 0; i < iterations; ++i, ++it) {
        size_t i_modded = (startpos + i) % v.size();
        if (v[i_modded] != *it) {
          RC_FAIL("v[" + std::to_string(i_modded) + "] = " +
                  std::to_string(v[i_modded]) + " == " + std::to_string(*it) +
                  " failed in forward.");
        }
      }

      // Test backward
      it = circular_begin(v, startpos);
      for (size_t i = 0; i < iterations; ++i, --it) {
        size_t i_modded = (v.size() - (i % v.size()) + startpos) % v.size();
        if (v[i_modded] != *it) {
          RC_FAIL("v[" + std::to_string(i_modded) + "] = " +
                  std::to_string(v[i_modded]) + " == " + std::to_string(*it) +
                  " failed in backward.");
        }
      }
    };

    REQUIRE(rc::check(
          "CircularIterator: Infinite circular iteration around vector", test));
  }

  //
  // ---------------------------------------------------------------
  //

  SECTION("Finite iteration around arbitrary vector with offset.") {
    auto test = [](std::vector<testtype> v) {
      size_t startpos = 0;
      if (v.size() > 0) {
        startpos = *gen::inRange<size_t>(0, v.size())
                          .as("position to start the circle");
      }

      // Construct the iterators:
      auto begin = circular_begin(v, startpos);
      auto end = circular_end(v, startpos);

      // Forward iteration:
      size_t n_iter = 0;
      for (; begin != end; ++begin) {
        size_t i_modded = (n_iter + startpos) % v.size();
        if (v[i_modded] != *begin) {
          RC_FAIL("v[" + std::to_string(i_modded) + "] = " +
                  std::to_string(v[i_modded]) + " == " +
                  std::to_string(*begin) + " failed in forward.");
        }
        ++n_iter;
      }
      RC_ASSERT(n_iter == v.size());

      // Backward iteration:
      n_iter = 0;
      begin = circular_begin(v, startpos);
      for (; begin != end; --begin) {
        size_t i_modded = (v.size() - n_iter + startpos) % v.size();
        if (v[i_modded] != *begin) {
          RC_FAIL("v[" + std::to_string(i_modded) + "] = " +
                  std::to_string(v[i_modded]) + " == " +
                  std::to_string(*begin) + " failed in backward.");
        }
        ++n_iter;
      }
      RC_ASSERT(n_iter == v.size());
    };

    REQUIRE(rc::check("CircularIterator: Iteration around vector with offset",
                      test));
  }

  //
  // ---------------------------------------------------------------
  //

  SECTION("Pushing elements back into circular buffer.") {
    auto test = [](std::vector<testtype> v) {
      RC_PRE(v.size() > 0u);

      CircularBuffer<testtype> buf{v.size() + 5};
      for (auto elem : v) {
        buf.push_back(elem);
      }

      // Access with iterators
      size_t i = 0;
      for (auto it = std::begin(buf); it != std::end(buf); ++it, ++i) {
        RC_ASSERT(*it == v[i]);
      }
      RC_ASSERT(i == v.size());

      // Access with const iterators
      i = 0;
      for (auto it = std::cbegin(buf); it != std::cend(buf); ++it, ++i) {
        RC_ASSERT(*it == v[i]);
      }
      RC_ASSERT(i == v.size());
    };

    REQUIRE(rc::check("Pushing elements into circular buffer", test));
  }

  //
  // ------------------------------------------------------------
  //

  SECTION("Random function test of circular buffer") {
    // Typedef the operations
    typedef PushBack<testtype> op_PushBack;
    typedef PushFront<testtype> op_PushFront;
    typedef Clear<testtype> op_Clear;
    typedef ChangeMaxSize<testtype> op_ChangeMaxSize;

    REQUIRE(rc::check(
          "Random function test of circular buffer with Clear",
          exectute_random_test<testtype, op_PushBack, op_PushFront, op_Clear>));
    REQUIRE(rc::check(
          "Random function test of circular buffer with ChangeMaxSize",
          exectute_random_test<testtype, op_PushBack, op_PushFront,
                               op_ChangeMaxSize>));
  }  // Random function test

}  // TEST_CASE
}  // namespace test
}  // namespace gscf
