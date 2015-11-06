//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_driver/Manager.hh
 * \author Thomas M>. Evans
 * \date   Wed Jun 18 11:21:16 2014
 * \brief  Manager class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_driver_Manager_hh
#define mc_driver_Manager_hh

#include <sstream>
#include <string>
#include <memory>

#include "comm/P_Stream.hh"
#include "mc/Fixed_Source_Solver.hh"
#include "mc/Keff_Solver.hh"
#include "mc/Global_RNG.hh"
#include "mc/Source_Transporter.hh"
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
    typedef profugus::Core                          Geom_t;
    typedef Problem_Builder::RCP_ParameterList      RCP_ParameterList;
    typedef Problem_Builder::SP_Physics             SP_Physics;
    typedef Problem_Builder::SP_Geometry            SP_Geometry;
    typedef profugus::Solver<Geom_t>                Solver_t;
    typedef std::shared_ptr<Solver_t>               SP_Solver;
    typedef profugus::Keff_Solver<Geom_t>           Keff_Solver_t;
    typedef std::shared_ptr<Keff_Solver_t>          SP_Keff_Solver;
    typedef profugus::Fixed_Source_Solver<Geom_t>   Fixed_Source_Solver_t;
    typedef std::shared_ptr<Fixed_Source_Solver_t>  SP_Fixed_Source_Solver;
    typedef Solver_t::Tallier_t                     Tallier_t;
    typedef Solver_t::SP_Tallier                    SP_Tallier;
    typedef profugus::Source_Transporter<Geom_t>    Transporter_t;
    typedef std::shared_ptr<Transporter_t>          SP_Transporter;
    typedef profugus::Global_RNG::RNG_Control_t     RNG_Control_t;
    typedef std::shared_ptr<RNG_Control_t>          SP_RNG_Control;
    typedef profugus::Fission_Source<Geom_t>        Fission_Source_t;
    typedef std::shared_ptr<Fission_Source_t>       SP_Fission_Source;

    // >>> DATA

    // Problem database.
    RCP_ParameterList d_db;

    // Geometry.
    SP_Geometry d_geometry;

    // Physics.
    SP_Physics d_physics;

    // Solvers.
    SP_Solver              d_solver;
    SP_Keff_Solver         d_keff_solver;
    SP_Fixed_Source_Solver d_fixed_solver;

    // Random number controller.
    SP_RNG_Control d_rnd_control;

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

    // Build an Anderson solver.
    template<class T>
    void build_anderson(SP_Transporter transporter, SP_Fission_Source source);

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

#endif // mc_driver_Manager_hh

//---------------------------------------------------------------------------//
//                 end of Manager.hh
//---------------------------------------------------------------------------//
