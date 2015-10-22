//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc_driver/omp.cc
 * \author Thomas M. Evans
 * \date   Wed Oct 21 10:43:31 2015
 * \brief  omp class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <string>
#include <iomanip>
#include <iostream>
#include <memory>
#include <cmath>

#include "comm/OMP.hh"
#include "comm/global.hh"
#include "Utils/rng/RNG_Control.hh"
#include "Utils/utils/Constants.hh"
#include "MC/geometry/Geometry.hh"
#include "MC/geometry/Mesh_Geometry.hh"
#include "MC/mc/Physics.hh"
#include "MC/mc/Particle.hh"
#include "MC/mc/Tallier.hh"
#include "MC/mc/Variance_Reduction.hh"
#include "MC/mc/VR_Roulette.hh"
#include "MC/mc/Domain_Transporter.hh"
#include "MC/mc/Uniform_Source.hh"
#include "MC/mc/Box_Shape.hh"

using Control    = profugus::RNG_Control;
using RNG        = Control::RNG_t;
using SP_Control = std::shared_ptr<Control>;

using Geometry    = profugus::Core;
using State       = Geometry::Geo_State_t;
using Core_t      = Geometry::Array_t;
using Lattice_t   = Core_t::Object_t;
using Pin_Cell_t  = Lattice_t::Object_t;
using SP_Core     = Geometry::SP_Array;
using SP_Lattice  = Core_t::SP_Object;
using SP_Pin_Cell = Lattice_t::SP_Object;
using SP_Geometry = std::shared_ptr<Geometry>;

using Mesh_Geometry = profugus::Mesh_Geometry;
using Mesh_State    = Mesh_Geometry::Geo_State_t;

using Physics_t  = profugus::Physics;
using SP_Physics = std::shared_ptr<Physics_t>;

using Bounds    = profugus::Group_Bounds;
using SP_Bounds = std::shared_ptr<Bounds>;

using XS     = Physics_t::XS_t;
using RCP_XS = Physics_t::RCP_XS;

using Tallier_t  = profugus::Tallier;
using SP_Tallier = std::shared_ptr<Tallier_t>;

using Var_Reduction_t  = profugus::Variance_Reduction;
using SP_Var_Reduction = std::shared_ptr<Var_Reduction_t>;

using Particle    = profugus::Particle;
using SP_Particle = std::shared_ptr<Particle>;

using Uniform_Source = profugus::Uniform_Source;
using Box            = profugus::Box_Shape;
using SP_Box         = std::shared_ptr<Box>;

using Space_Vector = def::Space_Vector;

const int N = 1000000;

using std::cout;
using std::endl;
using namespace profugus::constants;

extern RNG rng;
#pragma omp threadprivate(rng)
RNG rng;

//---------------------------------------------------------------------------//

Geometry build_geo()
{
    SP_Pin_Cell pin = std::make_shared<Pin_Cell_t>(0, 1.0, 1.0);
    SP_Lattice  lat = std::make_shared<Lattice_t>(10, 10, 10, 1);
    lat->assign_object(pin, 0);
    lat->complete(0, 0, 0);

    SP_Core core = std::make_shared<Core_t>(1, 1, 1, 1);
    core->assign_object(lat, 0);
    core->complete(0, 0, 0);

    Geometry geo(core);
    return geo;
}

//---------------------------------------------------------------------------//

Physics_t build_physics()
{
    RCP_XS xs(Teuchos::rcp(new XS()));
    xs->set(0, 1);

    XS::OneDArray tot(1, 1.0);
    XS::TwoDArray sct(1, 1, 0.4);
    XS::OneDArray bounds(def::Vec_Dbl{100, 0.1});

    xs->set_bounds(bounds);
    xs->add(0, XS::TOTAL, tot);
    xs->add(0, 0, sct);

    xs->complete();

    Physics_t::RCP_Std_DB db(
        Teuchos::rcp(new Physics_t::ParameterList_t()));

    Physics_t physics(db, xs);
    return physics;
}

//---------------------------------------------------------------------------//

double test_A()
{
    Control control(12523);

    // Build random numbers
#pragma omp parallel
    {
#pragma omp critical
        {
            rng = control.rng(profugus::thread_id());
        }
    }

    double sum = 0.0;

    double begin = profugus::thread_time();

#pragma omp parallel for reduction(+:sum)
    for (int n = 0; n < N; ++n)
    {
        sum += rng.ran();
    }

    double end = profugus::thread_time();

    cout << "Final result = " << sum / static_cast<double>(N)
         << endl;

    return end - begin;
}

//---------------------------------------------------------------------------//

double test_B()
{
    Control control(12523);

    // Build random numbers
#pragma omp parallel
    {
#pragma omp critical
        {
            rng = control.rng(profugus::thread_id());
        }
    }

    Geometry geo = build_geo();

    double begin = profugus::thread_time();

#pragma omp parallel for
    for (int n = 0; n < N; ++n)
    {
        // sample position
        Space_Vector pos;
        pos[0] = rng.ran() * 10.0;
        pos[1] = rng.ran() * 10.0;
        pos[2] = rng.ran() * 10.0;

        // sample direction
        Space_Vector dir;

        double  costheta = 1.0 - 2.0 * rng.ran();
        double  sintheta = std::sqrt(1.0 - costheta * costheta);
        double  phi      = two_pi * rng.ran();
        dir[0]           = sintheta * std::cos(phi);
        dir[1]           = sintheta * std::sin(phi);
        dir[2]           = costheta;

        // Initialize
        State state;
        geo.initialize(pos, dir, state);
    }

    double end = profugus::thread_time();

    return end - begin;
}

//---------------------------------------------------------------------------//

double test_C()
{
    Control control(12523);

    // Build random numbers
#pragma omp parallel
    {
#pragma omp critical
        {
            rng = control.rng(profugus::thread_id());
        }
    }

    def::Vec_Dbl r = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::shared_ptr<Mesh_Geometry> geo =
        std::make_shared<Mesh_Geometry>(r, r, r);

    double begin = profugus::thread_time();

#pragma omp parallel for
    for (int n = 0; n < N; ++n)
    {
        // sample position
        Space_Vector pos;
        pos[0] = rng.ran() * 10.0;
        pos[1] = rng.ran() * 10.0;
        pos[2] = rng.ran() * 10.0;

        // sample direction
        Space_Vector dir;

        double  costheta = 1.0 - 2.0 * rng.ran();
        double  sintheta = std::sqrt(1.0 - costheta * costheta);
        double  phi      = two_pi * rng.ran();
        dir[0]           = sintheta * std::cos(phi);
        dir[1]           = sintheta * std::sin(phi);
        dir[2]           = costheta;

        // Initialize
        Mesh_State state;
        geo->initialize(pos, dir, state);
    }

    double end = profugus::thread_time();

    return end - begin;
}

//---------------------------------------------------------------------------//

double test_D()
{
    Control control(12523);

    // Build random numbers
#pragma omp parallel
    {
#pragma omp critical
        {
            rng = control.rng(profugus::thread_id());
        }
    }

    SP_Geometry geo     = std::make_shared<Geometry>(build_geo());
    SP_Physics  physics = std::make_shared<Physics_t>(build_physics());
    physics->set_geometry(geo);

    SP_Tallier tallier = std::make_shared<Tallier_t>();
    tallier->set(geo, physics);
    tallier->build();

    Physics_t::RCP_Std_DB db(
        Teuchos::rcp(new Physics_t::ParameterList_t()));
    SP_Var_Reduction var = std::make_shared<profugus::VR_Roulette>(db);

    profugus::Domain_Transporter trans;
    trans.set(geo, physics);
    trans.set(var);
    trans.set(tallier);

    double begin = profugus::thread_time();

#pragma omp parallel
    {
        profugus::Domain_Transporter::Bank_t bank;

#pragma omp for
        for (int n = 0; n < N; ++n)
        {
            // sample position
            Space_Vector pos;
            pos[0] = 1.0 + rng.ran() * 9.0;
            pos[1] = 1.0 + rng.ran() * 9.0;
            pos[2] = 1.0 + rng.ran() * 9.0;

            // sample direction
            Space_Vector dir;

            double  costheta = 1.0 - 2.0 * rng.ran();
            double  sintheta = std::sqrt(1.0 - costheta * costheta);
            double  phi      = two_pi * rng.ran();
            dir[0]           = sintheta * std::cos(phi);
            dir[1]           = sintheta * std::sin(phi);
            dir[2]           = costheta;

            // Initialize particle
            SP_Particle p = std::make_shared<Particle>();
            p->set_wt(1.0);
            p->set_rng(rng);
            p->set_group(0);
            p->set_matid(0);
            p->live();

            // Initialize the geometry
            geo->initialize(pos, dir, p->geo_state());

            // Transport
            trans.transport(*p, bank);
        }
    }

    double end = profugus::thread_time();

    return end - begin;
}

//---------------------------------------------------------------------------//

double test_E()
{
    SP_Control  control = std::make_shared<Control>(12523);
    SP_Geometry geo     = std::make_shared<Geometry>(build_geo());
    SP_Physics  physics = std::make_shared<Physics_t>(build_physics());
    physics->set_geometry(geo);

    Physics_t::RCP_Std_DB db(
        Teuchos::rcp(new Physics_t::ParameterList_t()));
    SP_Var_Reduction var = std::make_shared<profugus::VR_Roulette>(db);

    db->set("Np", N);

    auto box = std::make_shared<Box>(1.0, 9.0, 1.0, 9.0, 1.0, 9.0);

    Uniform_Source source(db, geo, physics, control);
    source.build_source(box);
    CHECK(source.num_to_transport() == N);

    double begin = profugus::thread_time();

#pragma omp parallel
    {
        profugus::Domain_Transporter::Bank_t bank;

#pragma omp for
        for (int n = 0; n < N; ++n)
        {
            // Initialize particle
            SP_Particle p = source.get_particle();
        }
    }

    double end = profugus::thread_time();

    return end - begin;
}

//---------------------------------------------------------------------------//

double test_F()
{
    SP_Control  control = std::make_shared<Control>(12523);
    SP_Geometry geo     = std::make_shared<Geometry>(build_geo());
    SP_Physics  physics = std::make_shared<Physics_t>(build_physics());
    physics->set_geometry(geo);

    Physics_t::RCP_Std_DB db(
        Teuchos::rcp(new Physics_t::ParameterList_t()));
    SP_Var_Reduction var = std::make_shared<profugus::VR_Roulette>(db);

    db->set("Np", N);

    auto box = std::make_shared<Box>(1.0, 9.0, 1.0, 9.0, 1.0, 9.0);

    Uniform_Source source(db, geo, physics, control);
    source.build_source(box);
    CHECK(source.num_to_transport() == N);

    SP_Tallier tallier = std::make_shared<Tallier_t>();
    tallier->set(geo, physics);
    tallier->build();

    profugus::Domain_Transporter trans;
    trans.set(geo, physics);
    trans.set(var);
    trans.set(tallier);

    double begin = profugus::thread_time();

#pragma omp parallel
    {
        profugus::Domain_Transporter::Bank_t bank;

#pragma omp for
        for (int n = 0; n < N; ++n)
        {
            // Initialize particle
            SP_Particle p = source.get_particle();

            // Transport
            trans.transport(*p, bank);
        }
    }

    double end = profugus::thread_time();

    return end - begin;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    profugus::initialize(argc, argv);

    int num_threads = std::stoi(argv[1]);
    profugus::set_num_threads(num_threads);

    double begin = profugus::thread_time();

    double local = test_F();

    double end = profugus::thread_time();

    cout << profugus::num_available_threads() << " "
         << end - begin << " "
         << local
         << endl;

    profugus::finalize();
    return 0;
}

//---------------------------------------------------------------------------//
// end of MC/mc_driver/omp.cc
//---------------------------------------------------------------------------//
