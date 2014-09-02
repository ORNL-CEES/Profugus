//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Decomposition.cc
 * \author Thomas M. Evans
 * \date   Tue Aug 28 08:28:00 2007
 * \brief  Decomposition member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Decomposition.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for user-specified id ordering.
 *
 * \param local_to_global vector defined over the range [0,N_local) with
 * domain defined over [0,N_global) such that \c local_to_global[local_index]
 * \c = \c global_index
 */
Decomposition::Decomposition(const Vec_Int &local_to_global)
    :
#ifdef COMM_MPI
    d_comm(profugus::communicator)
#else
    d_comm()
#endif
    , d_map(Teuchos::rcp(new Map(-1, local_to_global.size(),
                                 &local_to_global[0], 0, d_comm)))
{
    REQUIRE(!d_map.is_null());

    INSIST(profugus::nodes() > 1, "Cannot construct with map on 1 pe.");

    ENSURE(d_map->NumMyElements() == local_to_global.size());
    ENSURE(d_comm.NumProc() == profugus::nodes());
    ENSURE(d_comm.MyPID() == profugus::node());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for user-specified id ordering on a pre-defined
 * communicator.
 *
 * \param local_to_global vector defined over the domain \c [0,N_local) with
 * range defined over \c [0,N_global) such that \c
 * local_to_global[local_index] \c = \c global_index and \c N_global is
 * defined relative to the number of processors in the communicator.
 *
 * \param comm communicator defining the processors that comprise the domain
 * \c [0,N_global)
 */
Decomposition::Decomposition(const Vec_Int        &local_to_global,
                             const Communicator_t &comm)
    :
#ifdef COMM_MPI
    d_comm(comm)
#else
    d_comm()
#endif
    , d_map(Teuchos::rcp(new Map(-1, local_to_global.size(),
                                 &local_to_global[0], 0, d_comm)))
{
    REQUIRE(!d_map.is_null());

    // set the local communicator
    profugus::set_internal_comm(comm);
    INSIST(profugus::nodes() > 1, "Cannot construct with map on 1 pe.");

    ENSURE(d_map->NumMyElements() == local_to_global.size());
    ENSURE(d_comm.NumProc() == profugus::nodes());
    ENSURE(d_comm.MyPID() == profugus::node());

    // set back to to Comm-world
    profugus::reset_internal_comm();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for global-id ordering (also useful for single-pe
 * problems).
 *
 * Pass the number of elements on each processor to the decomposition.  The
 * global-id is then defined:
 * \f[
   \mathrm{global_id} = \mathrm{local_id}_p + \mathrm{offset}_p\:,
 * \f]
 * where
 * \f[
   \mathrm{offset}_p = \sum_{i=0}^{p-1}\mathrm{num_elements}_i\:.
 * \f]
 *
 * \param num_elements number of local elements
 * \param local_map make the map local (i.e. fully replicated)
 */
Decomposition::Decomposition(int  num_elements,
                             bool local_map)
    :
#ifdef COMM_MPI
    d_comm(profugus::communicator)
#else
    d_comm()
#endif
{
    if( local_map )
        d_map = Teuchos::rcp(new LocalMap(num_elements, 0, d_comm));
    else
        d_map = Teuchos::rcp(new Map(-1, num_elements, 0, d_comm));

    ENSURE(profugus::nodes() == 1 || local_map ?
            d_map->NumGlobalElements() == num_elements :
            d_map->NumMyElements() < d_map->NumGlobalElements());
    ENSURE(d_map->NumMyElements() == num_elements);
    ENSURE(d_comm.NumProc() == profugus::nodes());
    ENSURE(d_comm.MyPID() == profugus::node());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for global-id ordering over a pre-defined communicator.
 *
 * Pass the number of elements on each processor to the decomposition.  The
 * difference here is that the communication space is defined by an input
 * communicator.  This allows clients to define space over a subset of total
 * processors that are determined by the communicator.  The global-id is then
 * defined:
 * \f[
   \mathrm{global_id} = \mathrm{local_id}_p + \mathrm{offset}_p\:,
 * \f]
 * where
 * \f[
   \mathrm{offset}_p = \sum_{i=0}^{p-1}\mathrm{num_elements}_i\:,
 * \f]
 * and the \c global_id is \b local to a communicator.
 *
 * \param num_elements number of local elements defined by the communicator
 * \param comm Communicator Object
 * \param local_map make the map local (i.e. fully replicated)
 */
Decomposition::Decomposition(int                   num_elements,
                             const Communicator_t &comm,
                             bool                  local_map)
    :
#ifdef COMM_MPI
    d_comm(comm)
#else
    d_comm()
#endif
{
    if( local_map )
        d_map = Teuchos::rcp(new LocalMap( num_elements, 0, d_comm));
    else
        d_map = Teuchos::rcp(new Map( -1, num_elements, 0, d_comm));

    REMEMBER(profugus::set_internal_comm(comm));
    ENSURE(profugus::nodes() == 1 || local_map ?
            d_map->NumGlobalElements() == num_elements :
            d_map->NumMyElements() < d_map->NumGlobalElements());
    ENSURE(d_map->NumMyElements() == num_elements);
    ENSURE(d_comm.NumProc() == profugus::nodes());
    ENSURE(d_comm.MyPID() == profugus::node());
    REMEMBER(profugus::reset_internal_comm());
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Decomposition.cc
//---------------------------------------------------------------------------//
