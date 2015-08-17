
//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaTypeTraits.hh
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_CudaTypeTraits_hh
#define mc_solver_CudaTypeTraits_hh

#include "CudaUtils.hh"

namespace alea{


template <typename T>
struct is_precomputed {
  static const bool value = false;
};
 
template <>
struct is_precomputed<OnTheFly> {
  static const bool value = false;
};
 
template <>
struct is_precomputed<Precomputed> {
  static const bool value = true;
};


}


#endif
