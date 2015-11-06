//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Tally.hh
 * \author Thomas M. Evans
 * \date   Wed May 14 15:10:09 2014
 * \brief  Tally base class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Tally_hh
#define mc_Tally_hh

#include <memory>
#include <string>

#include "Physics.hh"
#include "Definitions.hh"
#include "harness/DBC.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Tally
 * \brief Base class for tallies.
 */
//===========================================================================//

template <class Geometry>
class Tally
{
  public:
    //@{
    //! Typedefs.
    typedef Physics<Geometry>               Physics_t;
    typedef typename Physics_t::Particle_t  Particle_t;
    typedef std::shared_ptr<Physics_t>      SP_Physics;
    //@}

  protected:
    // >>> DATA

    // Physics.
    SP_Physics b_physics;

    // Tally name.
    std::string b_name;

    // Is this on during inactive cycles?
    bool b_inactive;

  public:
    //! Constructor.
    Tally(SP_Physics physics, bool inactive)
        : b_physics(physics)
        , b_inactive(inactive)
    {
        REQUIRE(b_physics);
    }

    // Destructor.
    virtual ~Tally() = 0;

    //! Set the tally name.
    void set_name(const std::string &name) { b_name = name; }

    //! Get the tally name.
    const std::string& name() const { return b_name; }

    //! Query if this tally is on during inactive cycles.
    bool inactive_cycle_tally() const { return b_inactive; }

    // >>> PUBLIC INTERFACE

    //! Accumulate first and second moments
    virtual void end_history() { /* * */ }

    //! Do post-processing on first and second moments
    virtual void finalize(double num_particles) { /* * */ }

    //! Begin active cycles in a kcode calculation (default no-op)
    virtual void begin_active_cycles() { /* * */ }

    //! Begin a new cycle in a kcode calculation (default no-op)
    virtual void begin_cycle() { /* * */ }

    //! End a cycle in a kcode calculation (default no-op)
    virtual void end_cycle(double num_particles) { /* * */ }

    //! Clear/re-initialize all tally values between solves
    virtual void reset() { /* * */ }
};

//---------------------------------------------------------------------------//
/*!
 * \class Source_Tally
 * \brief Defines source tally interfaces.
 */
template <class Geometry>
class Source_Tally : public Tally<Geometry>
{
    typedef Tally<Geometry>                 Base;
    typedef Physics<Geometry>               Physics_t;
    typedef typename Physics_t::Particle_t  Particle_t;
    typedef std::shared_ptr<Physics_t>      SP_Physics;

  public:
    // Constructor.
    Source_Tally(SP_Physics physics, bool inactive)
        : Base(physics, inactive)
    { /*...*/ }

    // Destructor.
    virtual ~Source_Tally() = 0;

    // >>> TALLY INTERFACE

    //! Tally events at particle birth.
    virtual void birth(const Particle_t &p) = 0;
};

//---------------------------------------------------------------------------//
/*!
 * \class Pathlength_Tally
 * \brief Defines source tally interfaces.
 */
template <class Geometry>
class Pathlength_Tally : public Tally<Geometry>
{
    typedef Tally<Geometry>                 Base;
    typedef Physics<Geometry>               Physics_t;
    typedef typename Physics_t::Particle_t  Particle_t;
    typedef std::shared_ptr<Physics_t>      SP_Physics;

  public:
    // Constructor.
    Pathlength_Tally(SP_Physics physics, bool inactive)
        : Base(physics, inactive)
    { /*...*/ }

    // Destructor.
    virtual ~Pathlength_Tally() = 0;

    // >>> TALLY INTERFACE

    //! Track particle and tally.
    virtual void accumulate(double step, const Particle_t &p) = 0;
};

//---------------------------------------------------------------------------//
/*!
 * \class Compound_Tally
 * \brief Tally that is multiple types (source and/or pathlength).
 */
template <class Geometry>
class Compound_Tally : public Tally<Geometry>
{
    typedef Tally<Geometry>                 Base;

  public:
    //@{
    //! Tally typedefs.aaa
    typedef Pathlength_Tally<Geometry>          Pathlength_Tally_t;
    typedef std::shared_ptr<Pathlength_Tally_t> SP_Pathlength_Tally;
    typedef Source_Tally<Geometry>              Source_Tally_t;
    typedef std::shared_ptr<Source_Tally_t>     SP_Source_Tally;
    typedef Physics<Geometry>                   Physics_t;
    typedef typename Physics_t::Particle_t      Particle_t;
    typedef std::shared_ptr<Physics_t>          SP_Physics;
    //@}

  protected:
    // >>> DATA

    // Tally components.
    SP_Pathlength_Tally b_pl_tally;
    SP_Source_Tally     b_src_tally;

  public:
    // Constructor.
    Compound_Tally(SP_Physics physics, bool inactive)
        : Base(physics, inactive)
    { /*...*/ }

    // Destructor.
    virtual ~Compound_Tally() = 0;

    // >>> TALLY INTERFACE

    //! Get the component pathlength tally.
    SP_Pathlength_Tally get_pl_tally() const { return b_pl_tally; }

    //! Get the component source tally.
    SP_Source_Tally get_src_tally() const { return b_src_tally; }
};

} // end namespace profugus

#endif // mc_Tally_hh

//---------------------------------------------------------------------------//
//                 end of Tally.hh
//---------------------------------------------------------------------------//
