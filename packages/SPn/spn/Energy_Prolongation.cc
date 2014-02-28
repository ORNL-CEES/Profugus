//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Energy_Prolongation.cc
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 12:35:44 2014
 * \brief  Energy_Prolongation member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "Epetra_Vector.h"

#include "harness/DBC.hh"
#include "Energy_Prolongation.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//

Energy_Prolongation::Energy_Prolongation(
        const Epetra_MultiVector &coarse_vec,
        const Epetra_MultiVector &fine_vec,
        const std::vector<int>   &steer_vec )
    : d_steer_vec(steer_vec)
    , d_coarse_map( coarse_vec.Map() )
    , d_fine_map( fine_vec.Map() )
{
    d_fine_groups = std::accumulate(steer_vec.begin(),steer_vec.end(),0);
    d_coarse_groups = steer_vec.size();
    Check( d_fine_groups > d_coarse_groups );

    // Determine energy-independent size of vector
    Check( fine_vec.MyLength()%d_fine_groups==0 );
    d_unks_per_grp   = fine_vec.MyLength() / d_fine_groups;
    Check( coarse_vec.MyLength()==d_unks_per_grp*d_coarse_groups );
}

//---------------------------------------------------------------------------//
// RESTRICTION OPERATOR
//---------------------------------------------------------------------------//

int Energy_Prolongation::Apply(
        const Epetra_MultiVector &coarse_vectors,
              Epetra_MultiVector &fine_vectors ) const
{
    Require( fine_vectors.MyLength()   == d_fine_groups*d_unks_per_grp );
    Require( coarse_vectors.MyLength() == d_coarse_groups*d_unks_per_grp );

    int num_vectors = fine_vectors.NumVectors();
    Check( coarse_vectors.NumVectors()==num_vectors );

    fine_vectors.PutScalar(0.0);

    // Process each vector in multivector
    int coarse_offset, fine_offset;
    for( int ivec=0; ivec<num_vectors; ++ivec )
    {
        // Access current vector
        const Epetra_Vector &coarse_vec = *(coarse_vectors(ivec));
        Epetra_Vector &fine_vec = *(fine_vectors(ivec));

        // Apply restriction to each component
        for( int i=0; i<d_unks_per_grp; ++i )
        {
            coarse_offset = i*d_coarse_groups;
            fine_offset   = i*d_fine_groups;
            int grp_ctr = 0;
            for( int icg=0; icg<d_coarse_groups; ++icg )
            {
                int fine_grps = d_steer_vec[icg];
                for( int ifg=grp_ctr; ifg<grp_ctr+fine_grps; ++ifg )
                {
                    fine_vec[fine_offset+ifg] = coarse_vec[coarse_offset+icg];
                }
                grp_ctr += fine_grps;
            }
        }
    }

    return 0;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Energy_Prolongation.cc
//---------------------------------------------------------------------------//
