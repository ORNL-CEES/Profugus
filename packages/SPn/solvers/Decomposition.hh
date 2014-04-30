//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Decomposition.hh
 * \author Thomas M. Evans
 * \date   Tue Aug 28 08:28:00 2007
 * \brief  Decomposition class definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_Decomposition_hh
#define solvers_Decomposition_hh

#include <SPn/config.h>
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Teuchos_RCP.hpp"

#include "harness/DBC.hh"
#include "comm/global.hh"
#include "utils/Definitions.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Decomposition
 * \brief Defines local-to-global mappings of a feature.
 *
 * For a given features, the Decomposition class defines local-to-global
 * mappings.  The feature could be cell-indices, node-indices, moment-indices,
 * etc.  The Decomposition class returns an Epetra_Map object that can then be
 * used to build Epetra objects for use with Trilinos solvers and services.
 *
 * The local-to-global mapping can be defined by the user, or simple global-id
 * ordering is applied.  See Decomposition() for more details.
 */
/*!
 * \example solvers/test/tstDecomposition.cc
 *
 * Test of Decomposition.
 */
//===========================================================================//

class Decomposition
{
  public:
    //! Integer-vector typedef.
    typedef def::Vec_Int             Vec_Int;
    typedef profugus::Communicator_t Communicator_t;

    //@{
    //! Epetra typedefs.
#ifdef COMM_MPI
    typedef Epetra_MpiComm    Comm;
#else
    typedef Epetra_SerialComm Comm;
#endif
    typedef Epetra_Map        Map;
    typedef Epetra_LocalMap   LocalMap;
    //@}
    typedef Teuchos::RCP<Map> RCP_Map;

  private:
    // >>> DATA

    // Communicator.
    Comm d_comm;

    // Map of local-to-global numbering.
    RCP_Map d_map;

  public:
    // Constructor with defined map.
    explicit Decomposition(const Vec_Int &local_to_global);

    // Defaults to global-id ordering.
    explicit Decomposition(int num_elements, bool local_map = false);

    // Constructor with defined map in a given communicator.
    Decomposition(const Vec_Int &local_to_global, const Communicator_t &comm);

    // Defaults to global-id ordering over a pre-defined communicator.
    Decomposition(int num_elements, const Communicator_t &comm,
                  bool local_map = false);

    // >>> ACCESSORS

    //! Get the communicator.
    const Comm& comm() const { return d_comm; }

    //! Get the map.
    const Map& map() const { return *d_map; }

    // Get the global index for a local index.
    inline int global(int local) const;

    //! Get the number of local elements on the node.
    int num_local() const { return d_map->NumMyElements(); }

    //! Get the number of global elements in the problem.
    int num_global() const { return d_map->NumGlobalElements(); }
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Get a global index for a given local index.
 */
int Decomposition::global(int local) const
{
    Require (local >= 0);
    Require (local < d_map->NumMyElements());

    return d_map->MyGlobalElements()[local];
}

} // end namespace profugus

#endif // solvers_Decomposition_hh

//---------------------------------------------------------------------------//
//              end of solvers/Decomposition.hh
//---------------------------------------------------------------------------//
