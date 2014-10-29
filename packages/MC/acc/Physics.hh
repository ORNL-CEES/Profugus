//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   acc/Physics.hh
 * \author Seth R Johnson
 * \date   Wed Oct 29 10:32:32 2014
 * \brief  Physics class declaration.
 * \note   Copyright (c) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef acc_Physics_hh
#define acc_Physics_hh


namespace profugus
{
class XS;
}

namespace acc
{

//===========================================================================//
/*!
 * \class Physics
 * \brief Flattend GPU physics class
 *
 * Long description or discussion goes here.
 *
 * \sa Physics.cc for detailed descriptions.
 *
 * \par Code Sample:
 * \code
 *     cout << "Hello, world." << endl;
 * \endcode
 */
/*!
 * \example Physics/test/tstPhysics.cc
 *
 * Test of Physics.
 */
//===========================================================================//

class Physics
{
  private:
    // >>> Device data

    // Total[matid][group]
    double *d_total;

    // Nusigf[matid][group]
    double *d_nusigf;

    // Scatter[matid][exiting][incident]
    double *d_scatter;

    // Number of materials/groups
    int d_num_mats;
    int d_num_groups;

    // >>> CPU

    // Memory storage on CPU
    std::vector<double> dv_total, dv_nusigf, dv_scatter;

  public:
    // Construct with number of matids, number of groups
    explicit Physics(const profugus::XS& xsdb);

  public:
    //! Number of elements in total, nusigf
    int num_vector_elements() const { return d_num_mat * d_num_groups; }

    //! Number of elements in the scattering matrix
    int num_matrix_elements() const { return d_num_mat * d_num_groups * d_num_groups; }

    //! Index into vector data
    int vector_index(int mat, int group) const
    {
        return mat * num_groups + group;
    }

    //! Index into matrix data
    int matrix_index(int mat, int out_group, int in_group) const
    {
        return (mat * num_groups + out_group) * num_groups + in_group;
    }

  private:
    // Copy data to GPU.
    void complete();
};

//---------------------------------------------------------------------------//
} // end namespace acc

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
// #include "Physics.i.hh"
//---------------------------------------------------------------------------//
#endif // acc_Physics_hh

//---------------------------------------------------------------------------//
// end of acc/Physics.hh
//---------------------------------------------------------------------------//
