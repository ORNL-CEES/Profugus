//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Device_Vector_device.pt.cc
 * \author Seth R Johnson
 * \date   Thu Aug 01 11:33:12 2013
 * \brief  Device_Vector explicit instantiation.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Device_Vector_device.t.hh"
#include "Definitions.hh"

namespace cuda
{
typedef arch::Device Arch_t;

template class Device_Vector<Arch_t,float>;
template class Device_Vector<Arch_t,double>;
template class Device_Vector<Arch_t,int>;
template class Device_Vector<Arch_t,unsigned int>;
template class Device_Vector<Arch_t,Space_Vector>;

template void device_to_host(
        const Device_Vector<Arch_t,float>&,
        profugus::View_Field<float>);
template void device_to_host(
        const Device_Vector<Arch_t,double>&,
        profugus::View_Field<double>);
template void device_to_host(
        const Device_Vector<Arch_t,int>&,
        profugus::View_Field<int>);
template void device_to_host(
        const Device_Vector<Arch_t,unsigned int>&,
        profugus::View_Field<unsigned int>);
template void device_to_host(
        const Device_Vector<Arch_t,Space_Vector>&,
        profugus::View_Field<Space_Vector>);

template void device_to_host(
        const Device_Vector<Arch_t,float>&,
        Host_Vector<float>&);
template void device_to_host(
        const Device_Vector<Arch_t,double>&,
        Host_Vector<double>&);
template void device_to_host(
        const Device_Vector<Arch_t,int>&,
        Host_Vector<int>&);
template void device_to_host(
        const Device_Vector<Arch_t,unsigned int>&,
        Host_Vector<unsigned int>&);
template void device_to_host(
        const Device_Vector<Arch_t,Space_Vector>&,
        Host_Vector<Space_Vector>&);

} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of Device_Vector_device.pt.cc
//---------------------------------------------------------------------------//
