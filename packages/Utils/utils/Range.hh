//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/utils/Range.hh
 * \author Konrad Rudolph and Seth R Johnson
 * \date   Sun Sep 20 09:55:26 2015
 * \brief  Range class declaration.
 * \note   Copyright (c) 2012 Konrad Rudolph, Licensed under Apache
 *
 * https://github.com/klmr/cpp11-range
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Range_hh
#define Utils_utils_Range_hh

#include <iterator>
#include <type_traits>

namespace profugus
{

//===========================================================================//
/*!
 * \fn range
 * \tparam T Value type to iterate over
 * \brief Get iterators over a range of values, or a semi-infinite range.
 *
 * \par Code Sample:
 * \code

    for (auto i : range(1, 5))
        cout << i << "\n";

    // Range of [0, 10)
    for (auto u : range(10u))
        cout << u << "\n";

    for (auto c : range('a', 'd'))
        cout << c << "\n";

    for (auto i : count(100).step(-3))
        if (i < 90) break;
        else        cout << i << "\n";

 * \endcode
 */
/*!
 * \example utils/test/tstRange.cc
 *
 * Test of range.
 */
//===========================================================================//

namespace detail
{
//===========================================================================//
/*!
 * \struct range_iter
 */
template<typename T>
struct range_iter : std::iterator<std::input_iterator_tag, T>
{
    range_iter(T current) : current(current) { }

    T operator *() const { return current; }

    T const* operator ->() const { return &current; }

    range_iter& operator ++()
    {
        ++current;
        return *this;
    }

    range_iter operator ++(int)
    {
        auto copy = *this;
        ++*this;
        return copy;
    }

    bool operator ==(range_iter const& other) const
    {
        return current == other.current;
    }

    bool operator !=(range_iter const& other) const
    {
        return !(*this == other);
    }

  protected:
    T current;
};

//===========================================================================//
/*!
 * \struct inf_range_iter
 *
 * Iterator that never finishes.
 */
template<typename T>
struct inf_range_iter : range_iter<T>
{
    inf_range_iter(T current = T()) : range_iter<T>(current) { }

    bool operator ==(inf_range_iter const&) const { return false; }
    bool operator !=(inf_range_iter const&) const { return true; }
};

//===========================================================================//
/*!
 * \struct step_range_iter
 */
template<typename T>
struct step_range_iter : range_iter<T>
{
    step_range_iter(T current, T step)
        : range_iter<T>(current), step(step) { }

    using range_iter<T>::current;

    step_range_iter& operator++()
    {
        current += step;
        return *this;
    }

    step_range_iter operator++(int)
    {
        auto copy = *this;
        ++*this;
        return copy;
    }

    // Loses commutativity. Iterator-based ranges are simply broken. :-(
    bool operator==(step_range_iter const& other) const
    {
        return step > 0 ? current >= other.current
            : current < other.current;
    }

    bool operator!=(step_range_iter const& other) const
    {
        return !(*this == other);
    }

  private:
    T step;
};

//===========================================================================//
/*!
 * \struct inf_step_range_iter
 *
 * Iterator that never finishes.
 */
template<typename T>
struct inf_step_range_iter : step_range_iter<T>
{
    inf_step_range_iter(T current = T(), T step = T())
        : step_range_iter<T>(current, step) { }

    bool operator ==(inf_step_range_iter const&) const { return false; }
    bool operator !=(inf_step_range_iter const&) const { return true; }
};

} // end profugus::detail

//===========================================================================//
/*!
 * \struct step_range_proxy
 *
 * Proxy container for iterating over a range of integral values with a step
 * between their values.
 */
template<typename T>
struct step_range_proxy
{
    using Iter_t = detail::step_range_iter<T>;

    step_range_proxy(T begin, T end, T step)
        : begin_(begin, step), end_(end, step) { }

    Iter_t begin() const { return begin_; }

    Iter_t end() const { return end_; }

  private:
    Iter_t begin_;
    Iter_t end_;
};

//===========================================================================//
/*!
 * \struct inf_step_range_proxy
 *
 * Proxy iterator for iterating over a range of integral values with a step
 * between their values.
 */
template<typename T>
struct inf_step_range_proxy
{
    using Iter_t = detail::inf_step_range_iter<T>;

    inf_step_range_proxy(T begin, T step) : begin_(begin, step) { }

    Iter_t begin() const { return begin_; }

    Iter_t end() const { return Iter_t(); }

  private:
    Iter_t begin_;
};

//===========================================================================//
/*!
 * \struct range_proxy
 *
 * Proxy iterator for iterating over a range of integral values.
 */
template <typename T>
struct range_proxy
{
    using Iter_t = detail::range_iter<T>;

    range_proxy(T begin, T end) : begin_(begin), end_(end) { }

    step_range_proxy<T> step(T step)
    {
        return {*begin_, *end_, step};
    }

    Iter_t begin() const { return begin_; }

    Iter_t end() const { return end_; }

  private:
    Iter_t begin_;
    Iter_t end_;
};

//===========================================================================//
/*!
 * \struct infinite_range_proxy
 *
 * Proxy iterator for iterating over a range of integral values.
 */
template <typename T>
struct infinite_range_proxy
{
    using Iter_t = detail::inf_range_iter<T>;

    infinite_range_proxy(T begin) : begin_(begin) { }

    inf_step_range_proxy<T> step(T step)
    {
        return {*begin_, step};
    }

    Iter_t begin() const { return begin_; }

    Iter_t end() const { return Iter_t(); }

  private:
    Iter_t begin_;
};

//===========================================================================//
/*!
 * \fn range
 *
 * Range over fixed beginning and end values.
 */
template <typename T>
range_proxy<T> range(T begin, T end)
{
    return {begin, end};
}

//===========================================================================//
/*!
 * \fn range
 *
 * Range with the default start value (0 for numeric types)
 */
template <typename T>
range_proxy<T> range(T end)
{
    return {T(), end};
}

//===========================================================================//
/*!
 * \fn count
 *
 * Count upward from zero.
 */
template <typename T>
infinite_range_proxy<T> count()
{
    return {T()};
}

//===========================================================================//
/*!
 * \fn count
 *
 * Count upward from a value.
 */
template <typename T>
infinite_range_proxy<T> count(T begin)
{
    return {begin};
}

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
#endif // Utils_utils_Range_hh

//---------------------------------------------------------------------------//
// end of Utils/utils/Range.hh
//---------------------------------------------------------------------------//
