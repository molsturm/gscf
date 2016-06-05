#pragma once
#include <iterator>
#include <linalgwrap/Exceptions.hh>
#include <list>
#include <memory>

// TODO This should go into the krikra library!
namespace gscf {

/** \name Class which allows to iterate between a range given by a start and end
 * iterator in circular fashion.
 *
 * This means that it may start at any point between m_first and m_last and it
 * iterates once through each element in this range, but flipping around from
 * m_last to m_first on calling operator++ */
template <typename Iterator>
class CircularIterator
      : public std::iterator<
              std::bidirectional_iterator_tag, typename Iterator::value_type,
              typename Iterator::difference_type, typename Iterator::pointer,
              typename Iterator::reference> {
public:
  typedef Iterator iterator_type;
  typedef typename iterator_type::value_type value_type;
  typedef typename iterator_type::reference reference;
  typedef typename iterator_type::pointer pointer;

  static_assert(std::is_same<typename iterator_type::iterator_category,
                             std::bidirectional_iterator_tag>::value ||
                      std::is_same<typename iterator_type::iterator_category,
                                   std::random_access_iterator_tag>::value,
                "The Iterator needs to be at least a bidirectional iterator");

  //@{
  /** Defaults for the big five */
  ~CircularIterator() = default;
  CircularIterator() = default;
  CircularIterator(CircularIterator&&) = default;
  CircularIterator(const CircularIterator&) = default;
  CircularIterator& operator=(CircularIterator&&) = default;
  CircularIterator& operator=(const CircularIterator&) = default;
  //@}

  /** \name Sensible construction yielding a valid CircularIterator
   *
   * Iterate on the range [first,end) in an infinite circle,
   * starting on position start.
   *
   * \param first  The first element of the range (should point to an actual
   * element!)
   * \param end    Iterator to past-the-last element (never touched)
   * \param start  A valid element within the range to start with.
   * \param begin_iterator Is this a begin iterator returned from a std::begin,
   * .begin() or similar function?
   *               This flag is needed to ensure the expected behaviour for
   * ranged-based iterators or iterators in for loops. The problem is that for
   * this circular iterator, the positon and range of the begin and end iterator
   * are all exactly identical.
   * (Since the range is a circle, the elmenet-past-the end is the beginning!)
   * */
  CircularIterator(iterator_type first, iterator_type end, iterator_type start,
                   bool begin_iterator = false);

  //
  // Increment and decrement
  //

  /** Prefix increment to the next element */
  CircularIterator& operator++();

  /** Postfix increment to the next element */
  CircularIterator operator++(int);

  /** Prefix increment to the next element */
  CircularIterator& operator--();

  /** Postfix increment to the next element */
  CircularIterator operator--(int);

  //
  // Element access
  //
  /** Return the value of the element we point to. */
  reference operator*() const;

  /** Access the members of the element we point to. */
  pointer operator->() const;

  //
  // Comparison
  //
  /** \brief check if two iterators are equal
   *
   * CircularIterators are equal if they iterate over the same range and point
   * to the same element and the begin_iterator flag is identical in both
   * objects.
   */
  bool operator==(const CircularIterator& other) const;

  /** \brief Check whether two iterators are unequal
   *
   * CircularIterators are unequal if they iterate over a different range or
   * point to a different element.
   */
  bool operator!=(const CircularIterator& other) const;

  //
  // Access range and position
  //
  /** Get the current iteration range */
  std::pair<iterator_type, iterator_type> iteration_range() const;

  /** Get the iterator to the current position */
  iterator_type position() const;

private:
  bool m_begin_iterator;  //< Is this a begin iterator
  iterator_type m_first;  //< First element of the range (inclusive)
  iterator_type m_pos;    //< Current position
  iterator_type m_end;    //< Last element of the range (exclusive)
};

/** \brief Convenience function to make a begin CircularIterator for a range,
 * specifying
 *  the start element as an iterator.
 *
 *  The iterator will iterate the range [begin:end), but in a circle and
 * starting
 *  from the element start.
 *  */
template <typename Iterator>
CircularIterator<Iterator> circular_begin(Iterator begin, Iterator end,
                                          Iterator start) {
  return CircularIterator<Iterator>(begin, end, start, true);
}

/** \brief Convenience function to make an end CircularIterator for a range,
 * specifying
 *  the start element as an iterator.
 *
 * The equivalent function to get a start iterator is circular_begin.
 * Use both function in pairs and supply them with the same arguments.
 *  */
template <typename Iterator>
CircularIterator<Iterator> circular_end(Iterator begin, Iterator end,
                                        Iterator start) {
  return CircularIterator<Iterator>(begin, end, start, false);
}

/** \brief Convenience function to make a begin CircularIterator for a
 * container, specifying the start element.
 *
 * The iterator will iterate the full container, but starting from the ith
 * element and looping around in circles.
 *  */
template <typename Container>
auto circular_begin(Container& c, typename Container::size_type i)
      -> CircularIterator<decltype(std::begin(c))> {
  return circular_begin(std::begin(c), std::end(c),
                        std::next(std::begin(c), i));
}

/** \brief Convenience function to make an end CircularIterator for a container,
 * specifying the start element.
 *
 * The equivalent function to get a start iterator is circular_begin.
 * Use both function in pairs and supply them with the same arguments.
 *  */
template <typename Container>
auto circular_end(Container& c, typename Container::size_type i)
      -> CircularIterator<decltype(std::begin(c))> {
  return circular_end(std::begin(c), std::end(c), std::next(std::begin(c), i));
}

/** \brief Convenience function to make a begin CircularIterator for a
 * container, specifying the start element.
 *
 * The iterator will iterate the full container, but starting from the ith
 * element and looping around in circles.
 *  */
template <typename Container>
auto circular_begin(Container& c, decltype(std::begin(c)) start)
      -> CircularIterator<decltype(std::begin(c))> {
  return circular_begin(std::begin(c), std::end(c), start);
}

/** \brief Convenience function to make an end CircularIterator for a container,
 * specifying the start element.
 *
 * The equivalent function to get a start iterator is circular_begin.
 * Use both function in pairs and supply them with the same arguments.
 *  */
template <typename Container>
auto circular_end(Container& c, decltype(std::begin(c)) start)
      -> CircularIterator<decltype(std::begin(c))> {
  return circular_end(std::begin(c), std::end(c), start);
}

//
// ----------------------------------------------------
//

template <typename T>
class CircularBuffer {
public:
  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;

  typedef std::list<T> container_type;
  typedef typename container_type::size_type size_type;
  typedef CircularIterator<typename container_type::iterator> iterator;
  typedef CircularIterator<typename container_type::const_iterator>
        const_iterator;

  /** \name Constructor
   *
   * \param max_size  The maximal size of the buffer
   * */
  CircularBuffer(size_type max_size);

  /** \name Constructor
   *
   * \param max_size  The maximal size of the buffer
   * \param il  initial list of elements for the buffer
   *            Assumes size of il to be no greater than
   *            max_size.
   **/
  CircularBuffer(size_type max_size, std::initializer_list<T> il);

  /* \name Modifiers
   */
  ///@{
  /** \name Push an element before the first element,
   *  possibly overwriting the current last element of the
   *  circular buffer if max_size() has been reached.
   */
  void push_front(value_type val);

  /** \name Push an element after the last element,
   *  possibly overwriting the current first element of the
   *  circular buffer if max_size() has been reached.
   */
  void push_back(value_type val);
  //@}

  /* \name Element access
   */
  ///@{
  //@{
  /* \name Access the first element */
  reference front();
  const_reference front() const;
  //@}

  //@{
  /* \name Access the last element */
  reference back();
  const_reference back() const;
  //@}
  ///@}

  /* \name Iterators
   */
  ///@{
  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;
  const_iterator cbegin() const;
  const_iterator cend() const;
  ///@}

  /** \name Discard all elements of the buffer
   *
   * The size is zero, but max_size is unaltered.
   * */
  void clear();
  ///@}

  /* \name Capacity
   */
  ///@{
  /** Return the maximal size of the buffer */
  size_type max_size() const;

  /** \brief Change the maximal size.
   *
   * If the max_size is increased, space for
   * more values is added at the back.
   * If max_size is decreased, leftover elements
   * at the back are deleted.
   */
  void max_size(size_type msize);

  /** Return the actual size of the buffer */
  size_type size() const;

  /** Test whether container is empty */
  bool empty() const;
  ///@}

private:
  //! Storage for the buffer
  container_type m_storage;

  //! Maximal size
  size_type m_max_size;

  /* Iterator to the virtual first element of the buffer
   *
   * Note that this is also the virtual past-the end element of this
   * (cyclic) buffer
   * */
  iterator m_first;
};

//
// ----------------------------------------------------------
//

template <typename Iterator>
CircularIterator<Iterator>::CircularIterator(iterator_type first,
                                             iterator_type end,
                                             iterator_type start,
                                             bool begin_iterator)
      : m_begin_iterator{begin_iterator},
        m_first{first},
        m_pos{start},
        m_end{end} {}

template <typename Iterator>
CircularIterator<Iterator>& CircularIterator<Iterator>::operator++() {
  assert_dbg(m_first != m_end,
             linalgwrap::ExcInvalidState(
                   "Cannot increment CircularIterator over empty range"));

  ++m_pos;
  if (m_pos == m_end) {
    // We are past the end: wrap around:
    m_pos = m_first;
  }

  m_begin_iterator = false;
  return *this;
}

template <typename Iterator>
CircularIterator<Iterator> CircularIterator<Iterator>::operator++(int) {
  CircularIterator copy{*this};
  ++(*this);
  return copy;
}

template <typename Iterator>
CircularIterator<Iterator>& CircularIterator<Iterator>::operator--() {
  assert_dbg(m_first != m_end,
             linalgwrap::ExcInvalidState(
                   "Cannot decrement CircularIterator over empty range"));

  if (m_pos == m_first) {
    m_pos = std::prev(m_end);
  } else {
    --m_pos;
  }

  m_begin_iterator = false;
  return *this;
}

template <typename Iterator>
CircularIterator<Iterator> CircularIterator<Iterator>::operator--(int) {
  CircularIterator copy{*this};
  --(*this);
  return copy;
}

template <typename Iterator>
typename CircularIterator<Iterator>::reference CircularIterator<Iterator>::
operator*() const {
  assert_dbg(m_first != m_end,
             linalgwrap::ExcInvalidState(
                   "Cannot dereference CircularIterator over empty range"));

  return *m_pos;
}

template <typename Iterator>
typename CircularIterator<Iterator>::pointer CircularIterator<Iterator>::
operator->() const {
  assert_dbg(m_first != m_end,
             linalgwrap::ExcInvalidState(
                   "Cannot dereference CircularIterator over empty range"));
  return m_pos.operator->();
}

template <typename Iterator>
bool CircularIterator<Iterator>::operator==(
      const CircularIterator& other) const {
  if (m_first == m_end && other.m_first == other.m_end &&
      m_first == other.m_first) {
    // Ranges of both iterators are equivalent and empty.
    // So m_begin_iterator and m_pos play no role.
    return true;
  }

  return m_pos == other.m_pos && m_first == other.m_first &&
         m_end == other.m_end && m_begin_iterator == other.m_begin_iterator;
}

template <typename Iterator>
bool CircularIterator<Iterator>::operator!=(
      const CircularIterator& other) const {
  return !(*this == other);
}

template <typename Iterator>
std::pair<typename CircularIterator<Iterator>::iterator_type,
          typename CircularIterator<Iterator>::iterator_type>
CircularIterator<Iterator>::iteration_range() const {
  return std::make_pair(m_first, m_end);
}

template <typename Iterator>
typename CircularIterator<Iterator>::iterator_type
CircularIterator<Iterator>::position() const {
  assert_dbg(m_first != m_end,
             linalgwrap::ExcInvalidState(
                   "Cannot get position iterator of CircularIterator "
                   "over empty range"));
  return m_pos;
}

//
// ---------------------------------------------------------------
//

template <typename T>
CircularBuffer<T>::CircularBuffer(size_type max_size)
      : CircularBuffer{max_size, std::initializer_list<T>{}} {}

template <typename T>
CircularBuffer<T>::CircularBuffer(size_type max_size,
                                  std::initializer_list<T> il)
      : m_storage{il},
        m_max_size{max_size},
        m_first{circular_begin(m_storage, 0)} {
  assert_greater_equal(il.size(), max_size);
}

template <typename T>
void CircularBuffer<T>::push_front(value_type val) {
  assert_dbg(max_size() != 0, linalgwrap::ExcInvalidState("max_size is zero"));

  if (m_storage.size() == m_max_size) {
    // Maximal size reached:
    // Decrement pointer (circularly) and change value
    --m_first;
    *m_first = std::move(val);
    return;
  }

  // Else we actually add an element:
  size_t ssize = m_storage.size();
  if (ssize == 0) {
    // Just push front and update iterator
    m_storage.push_front(std::move(val));
    m_first = circular_begin(m_storage, 0);
  } else {
    // Push and update m_first circular iterator.
    auto first_iterator = m_storage.insert(m_first.position(), std::move(val));
    m_first = circular_begin(m_storage, first_iterator);
  }

  // Check that we inserted something
  assert_dbg(ssize + 1 == m_storage.size(), linalgwrap::ExcInternalError());
}

template <typename T>
void CircularBuffer<T>::push_back(value_type val) {
  assert_dbg(max_size() != 0, linalgwrap::ExcInvalidState("max_size is zero"));

  if (m_storage.size() == m_max_size) {
    // Maximal size reached:
    // Change value and increment first pointer (circularly)
    // Here we use the fact that the past-the-end iterator in this circular
    // buffer is equvalant to m_first.
    *m_first = std::move(val);
    ++m_first;
    return;
  }

  size_t ssize = m_storage.size();
  if (ssize == 0 || m_first.position() == std::begin(m_storage)) {
    // if we are at the actual front, insert at back position
    m_storage.push_back(std::move(val));
    m_first = circular_begin(m_storage, 0);
  } else {
    // Else insert before the first element, i.e. circularly at the back
    m_storage.insert(m_first.position(), std::move(val));
    // Update circular iterator
    m_first = circular_begin(m_storage, m_first.position());
  }

  // Check that we inserted something
  assert_dbg(ssize + 1 == m_storage.size(), linalgwrap::ExcInternalError());
}

template <typename T>
typename CircularBuffer<T>::reference CircularBuffer<T>::front() {
  return *m_first;
}

template <typename T>
typename CircularBuffer<T>::const_reference CircularBuffer<T>::front() const {
  return *m_first;
}

template <typename T>
typename CircularBuffer<T>::reference CircularBuffer<T>::back() {
  // Return the last actual element of the buffer, this is the one
  // previous to first in the circular sense
  return *std::prev(m_first);
}

template <typename T>
typename CircularBuffer<T>::const_reference CircularBuffer<T>::back() const {
  // Return the last actual element of the buffer, this is the one
  // previous to first in the circular sense
  return *std::prev(m_first);
}

template <typename T>
typename CircularBuffer<T>::iterator CircularBuffer<T>::begin() {
  // Return a circular iterator, which is a begin iterator.
  return circular_begin<typename container_type::iterator>(
        std::begin(m_storage), std::end(m_storage), m_first.position());
}

template <typename T>
typename CircularBuffer<T>::iterator CircularBuffer<T>::end() {
  // Return a circular iterator, which is a begin iterator.
  return circular_end<typename container_type::iterator>(
        std::begin(m_storage), std::end(m_storage), m_first.position());
}

template <typename T>
typename CircularBuffer<T>::const_iterator CircularBuffer<T>::begin() const {
  return cbegin();
}

template <typename T>
typename CircularBuffer<T>::const_iterator CircularBuffer<T>::end() const {
  return cend();
}

template <typename T>
typename CircularBuffer<T>::const_iterator CircularBuffer<T>::cbegin() const {
  return circular_begin<typename container_type::const_iterator>(
        std::cbegin(m_storage), std::cend(m_storage), m_first.position());
}

template <typename T>
typename CircularBuffer<T>::const_iterator CircularBuffer<T>::cend() const {
  return circular_end<typename container_type::const_iterator>(
        std::cbegin(m_storage), std::cend(m_storage), m_first.position());
}

template <typename T>
void CircularBuffer<T>::clear() {
  // Clear storage and reset iterator:
  m_storage.clear();
  m_first = circular_begin(m_storage, 0);
}

template <typename T>
typename CircularBuffer<T>::size_type CircularBuffer<T>::max_size() const {
  return m_max_size;
}

template <typename T>
void CircularBuffer<T>::max_size(size_type msize) {
  m_max_size = msize;

  if (msize == 0) {
    // Clear all elements of the buffer:
    m_storage.clear();
    m_first = circular_begin(m_storage, 0);
  } else if (msize < m_storage.size()) {
    // We have to bin some elements of the buffer
    typedef typename container_type::iterator cont_iter;

    // The element after the one to be deleted is
    // the element past the last one in the new array.
    // In the circular sense, this is m_first
    iterator end_remove_range = m_first;

    // The first element we want to delete is the one that does
    // not fit inside the msize and is located the furthest away from the start.
    // This is the same as first advanced by msize times.
    iterator begin_remove_range = std::next(m_first, msize);

    // Delete elements at the back of the container.
    //
    // The idea is to first find the end of the first deletetion operaton,
    // by checking whatever comes first: end_remove_range or
    // std::end(m_storage)
    cont_iter del_begin = begin_remove_range.position();
    cont_iter del_end = del_begin;
    for (; del_end != std::end(m_storage) &&
           del_end != end_remove_range.position();
         ++del_end) {
      // Update the begin_remove_range,
      // such that we are prepared for the next step
      ++begin_remove_range;
    }
    m_storage.erase(del_begin, del_end);

    // Delete remaining elements which sit at the front
    del_begin = begin_remove_range.position();
    del_end = end_remove_range.position();
    m_storage.erase(del_begin, del_end);

    // Update m_first:
    m_first = circular_begin(m_storage, m_first.position());
  }
}

template <typename T>
typename CircularBuffer<T>::size_type CircularBuffer<T>::size() const {
  return m_storage.size();
}

template <typename T>
bool CircularBuffer<T>::empty() const {
  return m_storage.empty();
}

}  // namespace gscf
