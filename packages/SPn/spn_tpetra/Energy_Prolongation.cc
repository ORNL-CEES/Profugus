//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn_tpetra/Energy_Prolongation.cc
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 12:35:44 2014
 * \brief  Energy_Prolongation member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "harness/DBC.hh"
#include "Energy_Prolongation.hh"

namespace profugus
{
namespace tpetra
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//

Energy_Prolongation::Energy_Prolongation(
        Teuchos::RCP<const MV>  coarse_vec,
        Teuchos::RCP<const MV>  fine_vec,
        const std::vector<int> &steer_vec )
    : d_steer_vec(steer_vec)
    , d_coarse_map( coarse_vec->getMap() )
    , d_fine_map( fine_vec->getMap() )
{
    d_fine_groups = std::accumulate(steer_vec.begin(),steer_vec.end(),0);
    d_coarse_groups = steer_vec.size();
    Check( d_fine_groups > d_coarse_groups );

    // Determine energy-independent size of vector
    Check( fine_vec->getLocalLength()%d_fine_groups==0 );
    d_unks_per_grp   = fine_vec->getLocalLength() / d_fine_groups;
    Check( coarse_vec->getLocalLength()==d_unks_per_grp*d_coarse_groups );
}

//---------------------------------------------------------------------------//
// RESTRICTION OPERATOR
//---------------------------------------------------------------------------//

void Energy_Prolongation::apply( const MV &coarse_vectors, MV &fine_vectors,
        Teuchos::ETransp mode, double alpha, double beta ) const

{
    Require( alpha == Teuchos::ScalarTraits<double>::one() );
    Require( beta == Teuchos::ScalarTraits<double>::zero() );
    Require( fine_vectors.getLocalLength()   == d_fine_groups*d_unks_per_grp );
    Require( coarse_vectors.getLocalLength() ==
             d_coarse_groups*d_unks_per_grp );

    int num_vectors = fine_vectors.getNumVectors();
    Check( coarse_vectors.getNumVectors()==num_vectors );

    fine_vectors.putScalar(0.0);

    // Process each vector in multivector
    int coarse_offset, fine_offset;
    for( int ivec=0; ivec<num_vectors; ++ivec )
    {
        // Access current vector data
        Teuchos::ArrayRCP<const double> coarse_data =
            coarse_vectors.getData(ivec);
        Teuchos::ArrayRCP<double> fine_data =
            fine_vectors.getDataNonConst(ivec);

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
                    fine_data[fine_offset+ifg] =
                        coarse_data[coarse_offset+icg];
                }
                grp_ctr += fine_grps;
            }
        }
    }
}

} // end namespace profugus
} // end namespace tpetra

//---------------------------------------------------------------------------//
//                 end of Energy_Prolongation.cc
//---------------------------------------------------------------------------//
