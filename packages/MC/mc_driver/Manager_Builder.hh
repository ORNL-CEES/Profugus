//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_driver/Manager_Builder.hh
 * \author Steven Hamilton
 * \date   Wed Nov 25 11:25:01 2015
 * \brief  Manager_Builder class definition.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_driver_Manager_Builder_hh
#define mc_driver_Manager_Builder_hh

#include <memory>

#include "Manager_Base.hh"

namespace mc
{

//===========================================================================//
/*!
 * \class Manager_Builder
 * \brief Factory class for building templated Manager.
 */
//===========================================================================//

class Manager_Builder
{
  public:

      typedef std::shared_ptr<Manager_Base> SP_Manager_Base;

      static SP_Manager_Base build(const std::string &xml_file);

    //! Output messages in a common format.
#define SCREEN_MSG(stream)                            \
    {                                                 \
        std::ostringstream m;                         \
        m << ">>> " << stream;                        \
        profugus::pcout << m.str() << profugus::endl; \
    }
};

} // end namespace mc

#endif // mc_driver_Manager_Builder_hh

//---------------------------------------------------------------------------//
//                 end of Manager_Builder.hh
//---------------------------------------------------------------------------//
