
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
struct is_onthefly {
  static const bool value = false;
};
 
template <>
struct is_onthefly<OnTheFly> {
  static const bool value = true;
};
 
template <>
struct is_onthefly<Precomputed> {
  static const bool value = false;
};


}


#endif
