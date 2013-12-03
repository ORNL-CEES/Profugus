//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Comm_Traits.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 15:01:59 2008
 * \brief  Comm_Traits class definition.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef comm_Comm_Traits_hh
#define comm_Comm_Traits_hh

namespace profugus
{

//===========================================================================//
/*!
 * \struct Comm_Traits
 *
 * This struct and its specializations are used to implement the type-safe
 * default message tags in Comm.  Any other type-determined property needed in
 * Comm would also go here.
 */
//===========================================================================//

template<class T>
struct Comm_Traits
{
};

//---------------------------------------------------------------------------//
// SPECIALIZATION OF INTRINSIC ELEMENTAL TYPES
//---------------------------------------------------------------------------//

template<>
struct Comm_Traits<char>
{
    static const int tag = 431;
};

template<>
struct Comm_Traits<unsigned char>
{
    static const int tag = 432;
};

template<>
struct Comm_Traits<short>
{
    static const int tag = 433;
};

template<>
struct Comm_Traits<unsigned short>
{
    static const int tag = 434;
};

template<>
struct Comm_Traits<int>
{
    static const int tag = 435;
};

template<>
struct Comm_Traits<unsigned int>
{
    static const int tag = 436;
};

template<>
struct Comm_Traits<long>
{
    static const int tag = 437;
};

template<>
struct Comm_Traits<unsigned long>
{
    static const int tag = 438;
};

template<>
struct Comm_Traits<float>
{
    static const int tag = 439;
};

template<>
struct Comm_Traits<double>
{
    static const int tag = 440;
};

template<>
struct Comm_Traits<long double>
{
    static const int tag = 441;
};

//---------------------------------------------------------------------------//
// SPECIALIZATION OF INTRINSIC POINTER TYPES
//---------------------------------------------------------------------------//

template<>
struct Comm_Traits<char *>
{
    static const int tag = 451;
};

template<>
struct Comm_Traits<unsigned char *>
{
    static const int tag = 452;
};

template<>
struct Comm_Traits<short *>
{
    static const int tag = 453;
};

template<>
struct Comm_Traits<unsigned short *>
{
    static const int tag = 454;
};

template<>
struct Comm_Traits<int *>
{
    static const int tag = 455;
};

template<>
struct Comm_Traits<unsigned int *>
{
    static const int tag = 456;
};

template<>
struct Comm_Traits<long *>
{
    static const int tag = 457;
};

template<>
struct Comm_Traits<unsigned long *>
{
    static const int tag = 458;
};

template<>
struct Comm_Traits<float *>
{
    static const int tag = 459;
};

template<>
struct Comm_Traits<double *>
{
    static const int tag = 460;
};

template<>
struct Comm_Traits<long double *>
{
    static const int tag = 461;
};

} // end namespace profugus

#endif // comm_Comm_Traits_hh

//---------------------------------------------------------------------------//
//              end of comm/Comm_Traits.hh
//---------------------------------------------------------------------------//
