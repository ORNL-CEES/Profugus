//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   driver/Manager.hh
 * \author Thomas M>. Evans
 * \date   Wed Jun 18 11:21:16 2014
 * \brief  Manager class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef driver_Manager_hh
#define driver_Manager_hh

#include <sstream>
#include <string>

#include "comm/P_Stream.hh"
#include "Problem_Builder.hh"

namespace mc
{

//===========================================================================//
/*!
 * \class Manager
 * \brief Manager class that drives the MC miniapp.
 */
//===========================================================================//

class Manager
{
  private:
    // Typedefs.
    typedef Problem_Builder::RCP_ParameterList RCP_ParameterList;
    typedef Problem_Builder::SP_Physics        SP_Physics;
    typedef Problem_Builder::SP_Geometry       SP_Geometry;

    // >>> DATA

    // Problem database.
    RCP_ParameterList d_db;

    // Geometry.
    SP_Geometry d_geometry;

    // Physics.
    SP_Physics d_physics;

  public:
    // Constructor.
    Manager();

    // Setup the problem.
    void setup(const std::string &xml_file);

    // Solve the problem.
    void solve();

    // Output.
    void output();

  private:
    // >>> IMPLEMENTATION

    // Processor.
    int d_node, d_nodes;

    // Problem name.
    std::string d_problem_name;

    //! Output messages in a common format.
#define SCREEN_MSG(stream)                            \
    {                                                 \
        std::ostringstream m;                         \
        m << ">>> " << stream;                        \
        profugus::pcout << m.str() << profugus::endl; \
    }

};

} // end namespace mc

#endif // driver_Manager_hh

//---------------------------------------------------------------------------//
//                 end of Manager.hh
//---------------------------------------------------------------------------//
