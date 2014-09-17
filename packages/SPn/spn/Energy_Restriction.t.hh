//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Energy_Restriction.t.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 12:35:37 2014
 * \brief  Energy_Restriction template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Energy_Restriction_t_hh
#define spn_Energy_Restriction_t_hh

#include <algorithm>

#include "harness/DBC.hh"
#include "Energy_Restriction.hh"
#include "VectorTraits.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//

template <class T>
Energy_Restriction<T>::Energy_Restriction(
        Teuchos::RCP<const MAP> fine_map,
        Teuchos::RCP<const MAP> coarse_map,
        const std::vector<int> &steer_vec )
    : OperatorAdapter<T>(fine_map,coarse_map)
    , d_steer_vec(steer_vec)
    , d_fine_map( fine_map )
    , d_coarse_map( coarse_map )
{
    d_fine_groups = std::accumulate(steer_vec.begin(),steer_vec.end(),0);
    d_coarse_groups = steer_vec.size();
    CHECK( d_fine_groups > d_coarse_groups );

    // Determine energy-independent size of vector
    CHECK( VectorTraits<T>::local_size(fine_map)%d_fine_groups==0 );
    d_unks_per_grp = VectorTraits<T>::local_size(fine_map) / d_fine_groups;
    CHECK( VectorTraits<T>::local_size(coarse_map) ==
           d_unks_per_grp*d_coarse_groups );
}

//---------------------------------------------------------------------------//
// RESTRICTION OPERATOR
//---------------------------------------------------------------------------//

template <class T>
void Energy_Restriction<T>::ApplyImpl( const MV &fine_vectors,
                                             MV &coarse_vectors ) const
{
    REQUIRE( VectorTraits<T>::local_length(Teuchos::rcpFromRef(fine_vectors))
             == d_fine_groups*d_unks_per_grp );
    REQUIRE( VectorTraits<T>::local_length(Teuchos::rcpFromRef(coarse_vectors))
            == d_coarse_groups*d_unks_per_grp );

    int num_vectors = MVT::GetNumberVecs(fine_vectors);
    CHECK( MVT::GetNumberVecs(coarse_vectors) ==num_vectors );

    MVT::MvInit(coarse_vectors,0.0);

    // Process each vector in multivector
    int coarse_offset, fine_offset;
    for( int ivec=0; ivec<num_vectors; ++ivec )
    {
        Teuchos::ArrayView<const double> fine_data =
            VectorTraits<T>::get_data(Teuchos::rcpFromRef(fine_vectors),ivec);
        Teuchos::ArrayView<double> coarse_data =
            VectorTraits<T>::get_data_nonconst(
                Teuchos::rcpFromRef(coarse_vectors),ivec);

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
                    coarse_data[coarse_offset+icg] += fine_data[fine_offset+ifg];
                }
                grp_ctr += fine_grps;
                coarse_data[coarse_offset+icg] /=
                    static_cast<double>(fine_grps);
            }
        }
    }
}

} // end namespace profugus

#endif // spn_Energy_Restriction_t_hh

//---------------------------------------------------------------------------//
//                 end of Energy_Restriction.t.hh
//---------------------------------------------------------------------------//
