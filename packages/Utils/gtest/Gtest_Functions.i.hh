//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/gtest/Gtest_Functions.i.hh
 * \author Seth R Johnson
 * \date   Tue Dec 09 17:55:09 2014
 * \brief  Gtest_Functions inline definitions
 * \note   Copyright (c) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_gtest_Gtest_Functions_i_hh
#define Utils_gtest_Gtest_Functions_i_hh

#include <iterator>

#include "Utils/harness/Soft_Equivalence.hh"
#include "gtest.h"

namespace profugus
{

//---------------------------------------------------------------------------//
// HELPER CLASSES
//---------------------------------------------------------------------------//

template<class T>
struct NiceTypeName
{
    static std::string name() { return "UNKNOWN"; }
    static void print(std::ostream& os, const T& value) { os << value; }
};

template<>
struct NiceTypeName<int>
{
    static std::string name() { return "int"; }
    static void print(std::ostream& os, int value) { os << value; }
};

template<>
struct NiceTypeName<float>
{
    static std::string name() { return "float"; }
    static void print(std::ostream& os, float value) { os << value; }
};

template<>
struct NiceTypeName<double>
{
    static std::string name() { return "double"; }
    static void print(std::ostream& os, double value)
    {
        os << std::setprecision(8) << value;
    }
};

template<>
struct NiceTypeName<unsigned int>
{
    static std::string name() { return "unsigned int"; }
    static void print(std::ostream& os, unsigned int value)
    {
        os << value << 'u';
    }
};

template<>
struct NiceTypeName<std::size_t>
{
    static std::string name() { return "std::size_t"; }
    static void print(std::ostream& os, std::size_t value)
    {
        os << value << "ul";
    }
};

template<>
struct NiceTypeName<std::string>
{
    static std::string name() { return "std::string"; }

    static void print(std::ostream& os, const std::string& value)
    {
        os << '"' << value << '"';
    }
};

//---------------------------------------------------------------------------//
template<class Container>
struct ArrayValueTraits
{
    typedef typename Container::size_type  size_type;
    typedef typename std::decay<typename Container::value_type>::type
        value_type;
};


template<typename T, std::size_t N>
struct ArrayValueTraits< T[N] >
{
    typedef std::size_t size_type;
    typedef T           value_type;
};

template<typename Expected_T, typename Actual_T>
struct Failed_Value
{
    typedef std::size_t size_type;
    typedef Expected_T  value_type_E;
    typedef Actual_T    value_type_A;

    size_type  index;
    value_type_E expected;
    value_type_A actual;
};

//---------------------------------------------------------------------------//
// TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Print expected values.
 */
template<class Container>
void print_expected(const Container& data, const char* label)
{
    typedef typename ArrayValueTraits<Container>::value_type value_type;
    using NiceType_t = NiceTypeName<value_type>;
    using std::cout;
    using std::endl;

    cout << "const " << NiceTypeName<value_type>::name()
         << " expected_" << label << "[] = {";
    auto iter = std::begin(data);
    auto end_iter = std::end(data);
    if (iter != end_iter)
    {
        NiceType_t::print(cout, *iter++);
    }
    while (iter != end_iter)
    {
        cout << ", ";
        NiceType_t::print(cout, *iter++);
    }
    cout << "};" << std::setprecision(6) << endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Vector comparison.
 */
template<class Container_E, class Container_A>
::testing::AssertionResult IsVecEq(
        const char* expected_expr,
        const char* actual_expr,
        const Container_E& expected,
        const Container_A& actual)
{
    typedef typename ArrayValueTraits<Container_E>::size_type  size_type;
    typedef typename ArrayValueTraits<Container_E>::value_type value_type_E;
    typedef typename ArrayValueTraits<Container_A>::value_type value_type_A;
    typedef Failed_Value<value_type_E, value_type_A>           Failed_Value_t;
    typedef std::vector<Failed_Value_t>                        Vec_Failed;

    size_type expected_size = std::end(expected) - std::begin(expected);
    size_type actual_size   = std::end(actual) - std::begin(actual);

    // First, check that the sizes are equal
    if (expected_size != actual_size)
    {
        ::testing::AssertionResult failmsg = ::testing::AssertionFailure();

        failmsg << " Size of: " << actual_expr << "\n"
                << "  Actual: " << actual_size << "\n"
                << "Expected: " << expected_expr << ".size()\n"
                << "Which is: " << expected_size << "\n";
        return failmsg;
    }

    // Keep track of what elements failed
    Vec_Failed failures;

    // Now loop through elements of the vectors and check for soft equivalence
    auto expected_iter = std::begin(expected);
    auto expected_end  = std::end(expected);
    auto actual_iter   = std::begin(actual);
    size_type i = 0;
    for (; expected_iter != expected_end; ++expected_iter, ++actual_iter, ++i)
    {
        if (*actual_iter != *expected_iter)
        {
            Failed_Value_t temp = {i, *expected_iter, *actual_iter};
            failures.push_back(temp);
        }
    }

    if (failures.empty())
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        using std::setw;
        using std::setprecision;

        ::testing::AssertionResult failmsg = ::testing::AssertionFailure();

        failmsg << "Values in: " << actual_expr << "\n"
                << " Expected: " << expected_expr << "\n"
                << failures.size() << " of "
                << expected_size << " elements differ:\n";

        // Only print the first 40 failures
        auto end_failures = failures.end();
        if (failures.size() > 40)
        {
            failmsg << "(Truncating to first 40 failed values)\n";
            end_failures = failures.begin() + 40;
        }

        // Calculate how many digits we need to space out
        unsigned int num_digits = calc_num_digits(failures.back().index);

        // Construct our own stringstream because google test ignores setw
        std::ostringstream failure_stream;

        // Try to use user-given expressions for headers, but fall back if the
        // column length is exceeded
        std::string e_expr(expected_expr);
        std::string a_expr(actual_expr);

        failure_stream
            << setw(num_digits) << "i" << " "
            << setw(16) << (e_expr.size() <= 16 ? e_expr : "EXPECTED") << " "
            << setw(16) << (a_expr.size() <= 16 ? a_expr : "ACTUAL") << "\n";

        // Loop through failed indices and print values
        for (auto it = failures.begin(); it != end_failures; ++it)
        {
            failure_stream
                << setw(num_digits) << it->index << " "
                << setw(16) << it->expected << " "
                << setw(16) << it->actual << "\n";
        }
        failmsg << failure_stream.str();

        return failmsg;
    }
}
//-------------------------------------------------------------------------//
/*!
 * \brief Custom vector comparison with soft equiavelence
 */
template<class Container_E, class Container_A>
::testing::AssertionResult IsVecSoftEquiv(
        const char* expected_expr,
        const char* actual_expr,
        const char* rel_expr,
        const char* abs_expr,
        const Container_E& expected,
        const Container_A& actual,
        double rel,
        double abs)
{
    typedef typename ArrayValueTraits<Container_E>::size_type  size_type;
    typedef typename ArrayValueTraits<Container_E>::value_type value_type_E;
    typedef typename ArrayValueTraits<Container_A>::value_type value_type_A;
    typedef std::numeric_limits<value_type_E>                  limits_t;
    typedef Failed_Value<value_type_E, value_type_A>           Failed_Value_t;
    typedef std::vector<Failed_Value_t>                        Vec_Failed;

    // Uncomment if you want to include type_traits
#if 0
    static_assert(std::is_convertible<
                        typename ArrayValueTraits<Container_A>::value_type,
                        value_type>,
            "Incompatible container classes given.");
#endif

    size_type expected_size = std::end(expected) - std::begin(expected);
    size_type actual_size   = std::end(actual) - std::begin(actual);

    // First, check that the sizes are equal
    if (expected_size != actual_size)
    {
        ::testing::AssertionResult failmsg = ::testing::AssertionFailure();

        failmsg << " Size of: " << actual_expr << "\n"
                << "  Actual: " << actual_size << "\n"
                << "Expected: " << expected_expr << ".size()\n"
                << "Which is: " << expected_size << "\n";
        return failmsg;
    }

    // Now loop through elements of the vectors and check for soft equivalence
    // Keep track of what elements failed
    Vec_Failed failures;

    // Now loop through elements of the vectors and check for soft equivalence
    auto expected_iter = std::begin(expected);
    auto expected_end  = std::end(expected);
    auto actual_iter   = std::begin(actual);
    size_type i = 0;
    for (; expected_iter != expected_end; ++expected_iter, ++actual_iter, ++i)
    {
        if (std::isinf(*expected_iter) && std::isinf(*actual_iter))
        {
            if (std::signbit(*expected_iter) == std::signbit(*actual_iter))
            {
                // Signs of infinities match
                continue;
            }
        }
        else if (profugus::soft_equiv(*actual_iter, *expected_iter, rel, abs))
        {
            continue;
        }

        // Failure
        Failed_Value_t temp = {i, *expected_iter, *actual_iter};
        failures.push_back(temp);
    }

    if (failures.empty())
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        using std::setw;
        using std::setprecision;

        ::testing::AssertionResult failmsg = ::testing::AssertionFailure();

        failmsg << "Values in: " << actual_expr << "\n"
                << " Expected: " << expected_expr << "\n"
                << failures.size() << " of "
                << expected_size << " elements differ more than "
                << "given tolerance " << rel_expr << ":\n";

        // Only print the first 40 failures
        auto end_failures = failures.end();
        if (failures.size() > 40)
        {
            failmsg << "(Truncating to first 40 failed values)\n";
            end_failures = failures.begin() + 40;
        }

        // Calculate how many digits we need to space out
        unsigned int num_digits = 0;
        {
            size_type temp = failures.back().index;
            do {
                temp /= 10;
                ++num_digits;
            } while (temp > 0);
        }

        // Construct our own stringstream because google test ignores setw
        std::ostringstream failure_stream;
        failure_stream << setprecision(limits_t::digits10);

        value_type_E error = -1;

        // Try to use user-given expressions for headers, but fall back if the
        // column length is exceeded
        std::string e_expr(expected_expr);
        std::string a_expr(actual_expr);

        failure_stream
            << setw(num_digits) << "i" << " "
            << setw(16) << (e_expr.size() <= 16 ? e_expr : "EXPECTED") << " "
            << setw(16) << (a_expr.size() <= 16 ? a_expr : "ACTUAL") << " "
            << setw(16) << "Difference" << "\n";

        // Loop through failed indices and print values
        for (auto it = failures.begin(); it != end_failures; ++it)
        {
            if (std::isinf(it->expected))
            {
                error = std::numeric_limits<value_type_E>::infinity();
            }
            else if (std::fabs(it->expected) > abs)
            {
                error = (it->actual - it->expected) / it->expected;
            }
            else
            {
                error = it->actual - it->expected;
            }

            failure_stream
                << setw(num_digits) << it->index << " "
                << setw(16) << it->expected << " "
                << setw(16) << it->actual << " "
                << setw(16) << error << "\n";
        }
        failmsg << failure_stream.str();

        return failmsg;
    }
}

//-------------------------------------------------------------------------//
/*!
 * \brief Custom vector comparison with default soft equiavelence
 *
 * This signature uses the default tolerance for the appropriate floating point
 * operations.
 */
template<class Container_E, class Container_A>
::testing::AssertionResult IsVecSoftEquiv(
        const char* expected_expr,
        const char* actual_expr,
        const Container_E& expected,
        const Container_A& actual)
{
    typedef typename ArrayValueTraits<Container_E>::value_type value_type_E;
    typedef typename ArrayValueTraits<Container_A>::value_type value_type_A;
    typedef typename detail::precision_type<value_type_E,value_type_A>::type
            precision_type;

    // Get default precision
    const precision_type rel
        = detail::softeq_traits<precision_type>::rel_prec();
    const precision_type abs
        = detail::softeq_traits<precision_type>::abs_thresh();

    return IsVecSoftEquiv(
            expected_expr, actual_expr,
            std::to_string(rel).c_str(), std::to_string(abs).c_str(),
            expected, actual, rel, abs);
}

//-------------------------------------------------------------------------//
/*!
 * \brief Custom vector comparison with default soft equiavelence
 *
 * This signature uses the default tolerance for the appropriate floating point
 * operations.
 */
template<class Container_E, class Container_A>
::testing::AssertionResult IsVecSoftEquiv(
        const char* expected_expr,
        const char* actual_expr,
        const char* rel_expr,
        const Container_E& expected,
        const Container_A& actual,
        const double rel)
{
    typedef typename ArrayValueTraits<Container_E>::value_type value_type_E;
    typedef typename ArrayValueTraits<Container_A>::value_type value_type_A;
    typedef typename detail::precision_type<value_type_E,value_type_A>::type
            precision_type;

    // Get default precision
    const precision_type abs
        = detail::softeq_traits<precision_type>::abs_thresh();

    return IsVecSoftEquiv(
            expected_expr, actual_expr,
            rel_expr, std::to_string(abs).c_str(),
            expected, actual, rel, abs);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Custom error mesages for relative error soft equiavelence
 */
template<class Value_E, class Value_A>
::testing::AssertionResult IsSoftEquiv(
        const char* expected_expr,
        const char* actual_expr,
        const char* rel_expr,
        const char* abs_expr,
        Value_E expected,
        Value_A actual,
        double rel,
        double abs)
{
    typedef typename detail::precision_type<Value_E,Value_A>::type
            precision_type;

    // Check for infinities
    if (std::isinf(expected) && std::isinf(actual))
    {
        if (std::signbit(expected) == std::signbit(actual))
        {
            // Signs of infinities match
            return ::testing::AssertionSuccess();
        }
    }

    // Normal numbers
    if (profugus::soft_equiv<Value_E, Value_A>(actual, expected, rel, abs))
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        ::testing::AssertionResult failure = ::testing::AssertionFailure();

        failure << "Value of: " << actual_expr << "\n"
                << "  Actual: " << actual << "\n"
                << "Expected: " << expected_expr << "\n"
                << "Which is: " << expected << "\n";

        if (std::fabs(static_cast<precision_type>(expected))
            < static_cast<precision_type>(abs))
        {
            // Avoid divide by zero errors
            failure << "(Absolute error " << actual - expected
                    << " exceeds tolerance " << abs_expr << ")";
        }
        else
        {
            failure << "(Relative error " << (actual - expected) / expected
                    << " exceeds tolerance " << rel_expr << ")";
        }

        return failure;
    }
}

//-------------------------------------------------------------------------//
/*!
 * \brief Soft equiavelence with default absolute error
 */
template<class Value_E, class Value_A>
::testing::AssertionResult IsSoftEquiv(
        const char* expected_expr,
        const char* actual_expr,
        const char* rel_expr,
        Value_E expected,
        Value_A actual,
        double rel)
{
    typedef typename detail::precision_type<Value_E,Value_A>::type
            precision_type;

    // Get default precision
    const precision_type abs
        = detail::softeq_traits<precision_type>::abs_thresh();

    return IsSoftEquiv(
            expected_expr, actual_expr,
            rel_expr, std::to_string(abs).c_str(),
            expected, actual, rel, abs);
}

//-------------------------------------------------------------------------//
/*!
 * \brief Soft equiavelence with default absolute error
 */
template<class Value_E, class Value_A>
::testing::AssertionResult IsSoftEquiv(
        const char* expected_expr,
        const char* actual_expr,
        Value_E expected,
        Value_A actual)
{
    typedef typename detail::precision_type<Value_E,Value_A>::type
            precision_type;

    // Get default precision
    const precision_type rel
        = detail::softeq_traits<precision_type>::rel_prec();
    const precision_type abs
        = detail::softeq_traits<precision_type>::abs_thresh();

    return IsSoftEquiv(
            expected_expr, actual_expr,
            std::to_string(rel).c_str(), std::to_string(abs).c_str(),
            expected, actual, rel, abs);
}

//---------------------------------------------------------------------------//
} // end namespace profugus

#endif // Utils_gtest_Gtest_Functions_i_hh

//---------------------------------------------------------------------------//
// end of gtest/Gtest_Functions.i.hh
//---------------------------------------------------------------------------//
