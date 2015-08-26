//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/comm/OMP.hh
 * \author Thomas M. Evans
 * \date   Thu Aug 20 10:18:57 2015
 * \brief  OMP class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_comm_OMP_hh
#define Utils_comm_OMP_hh

#ifdef _OPENMP
#include <omp.h>
#endif

namespace profugus
{

//---------------------------------------------------------------------------//
// OpenMP FUNCTIONS
//---------------------------------------------------------------------------//

inline bool multithreading_available()
{
#ifdef _OPENMP
    return true;
#endif
    return false;
}

//---------------------------------------------------------------------------//

inline void set_num_threads(int nt)
{
#ifdef _OPENMP
    omp_set_num_threads(nt);
#endif
}

//---------------------------------------------------------------------------//

inline int num_current_threads()
{
#ifdef _OPENMP
    return omp_get_num_threads();
#endif
    return 1;
}

//---------------------------------------------------------------------------//

inline int thread_id()
{
#ifdef _OPENMP
    return omp_get_thread_num();
#endif
    return 0;
}

//---------------------------------------------------------------------------//

inline int num_available_threads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#endif
    return 1;
}

//---------------------------------------------------------------------------//

inline void turn_off_dynamic_threading()
{
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
}

//---------------------------------------------------------------------------//

inline void turn_on_dynamic_threading()
{
#ifdef _OPENMP
    omp_set_dynamic(1);
#endif
}

//---------------------------------------------------------------------------//

inline bool is_threading_dynamic()
{
#ifdef _OPENMP
    return omp_get_dynamic();
#endif
    return false;
}

//---------------------------------------------------------------------------//

inline bool in_thread_parallel_region()
{
#ifdef _OPENMP
    return omp_in_parallel();
#endif
    return false;
}

} // end namespace profugus

#endif // Utils_comm_OMP_hh

//---------------------------------------------------------------------------//
// end of Utils/comm/OMP.hh
//---------------------------------------------------------------------------//
