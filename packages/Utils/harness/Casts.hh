//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/Casts.hh
 * \author Seth R Johnson
 * \date   Wed Sep 18 14:20:31 2013
 * \brief  Inline implementation for different casts
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef harness_Casts_hh
#define harness_Casts_hh

#include "DBC.hh"

namespace nemesis
{
//---------------------------------------------------------------------------//
// Only declare base class; don't define
template<typename Derived>
class smart_cast;

/*!
 * \brief Perform a valid cast on a pointer without a dynamic_cast
 *
 * This allows safely passing around base-class pointers between non-templated
 * code, then inside a template derived class extracting the corresponding
 * templated pointer/reference.
 *
 * Using a static_cast instead of a reinterpret_cast causes a compile-time
 * error if two objects don't share the same class hierarchy.
 *
 * \warning The base class must polymorphic (virtual functions)
 *
 * Based on code at
 * http://www.drdobbs.com/templates-for-efficient-dynamic-type-che/184403724?pgno=3
 */
template<typename Derived>
class smart_cast<Derived*>
{
  public:
    template<typename Base>
    smart_cast(Base* base_ptr)
      : d_result(static_cast<Derived*>(base_ptr))
    {
        Require(dynamic_cast<Derived*>(base_ptr) == d_result);
    }

    operator Derived*() const { return d_result; }

  private:
    Derived* d_result;
};

//! Specialization for references
template<typename Derived>
class smart_cast<Derived&>
{
  public:
    template<typename Base>
    smart_cast(Base& base_ref)
      : d_result(static_cast<Derived&>(base_ref))
    {
        Require(dynamic_cast<Derived*>(&base_ref) == &d_result);
    }

    operator Derived&() const { return d_result; }

  private:
    Derived& d_result;
};
//---------------------------------------------------------------------------//
} // end namespace nemesis

#endif // harness_Casts_hh

//---------------------------------------------------------------------------//
//                 end of Casts.hh
//---------------------------------------------------------------------------//
