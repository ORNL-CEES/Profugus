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

#include <vector>
#include "Particle.hh"

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

    // Scatter[matid][incident][exiting]
    double *d_outscatter_pdf;

    // Scattering ratio[matid][group]
    double *d_scatter_ratio;

    // Fissionable[matid]
    int *d_fissionable;

    // Number of materials/groups
    int d_num_mats;
    int d_num_groups;

    // >>> CPU

    // Memory storage on CPU
    std::vector<double> dv_total, dv_nusigf, dv_outscatter_pdf, dv_scatter_ratio;
    std::vector<int> dv_fissionable;

  public:
    // Construct with number of matids, number of groups
    explicit Physics(const profugus::XS& xsdb);

    // Clean up memory
    ~Physics();

  public:
    //! Total macro XS
#pragma acc routine seq
    double total(int matid, int group) const;

    //! Scattering ratio
#pragma acc routine seq
    double scattering_ratio(int matid, int group) const;

    //! Nu fission
#pragma acc routine seq
    double nusigf(int matid, int group) const;

    //! Collision
#pragma acc routine seq
    void collide(Particle& p) const;

    //! Is the material fissionable?
    bool is_fissionable(int mat) const
    {
        return d_fissionable[mat];
    }

  public:
    //! Number of elements in total, nusigf
    int num_vector_elements() const { return d_num_mats * d_num_groups; }

    //! Number of elements in the scattering matrix
    int num_matrix_elements() const { return d_num_mats * d_num_groups * d_num_groups; }

    //! Index into vector data
#pragma acc routine seq
    int vector_index(int mat, int group) const;

    //! Index into matrix data
#pragma acc routine seq
    int matrix_index(int mat, int out_group, int in_group) const;
  private:
    // Copy data to GPU.
    void complete();
};

//---------------------------------------------------------------------------//
} // end namespace acc
#endif // acc_Physics_hh

//---------------------------------------------------------------------------//
// end of acc/Physics.hh
//---------------------------------------------------------------------------//
