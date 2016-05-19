//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Bank.hh
 * \author Stuart Slattery
 * \brief  Bank class definition
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Bank_hh
#define cuda_mc_Bank_hh

#include <memory>
#include <vector>

#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Shared_Device_Ptr.hh"

#include "Definitions.hh"
#include "Particle.hh"
#include "Particle_Vector.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Bank
 * \brief Base class for storage for particles
 *
 * The bank is essentially a stack of particles. To reduce storage, we use a
 * delayed-copy approach.
 */
/*!
 * \example mc/test/tstBank.cc
 *
 * Test of Bank.
 */

//===========================================================================//
template <class Geometry>
class Bank
{
  public:
    //@{
    //! Useful typedefs.
    typedef Particle<Geometry>           Particle_t;
    typedef std::shared_ptr<Particle_t>  SP_Particle;
    typedef SP_Particle                  value_type;
    typedef size_t                       size_type;
    typedef SP_Particle&                 reference;
    typedef const SP_Particle&           const_reference;
    //@}

  private:
    // Container type for particles
    typedef std::vector<SP_Particle> Stack_Particle;

    // Container type for number of copies per particle
    typedef std::vector<size_type> Stack_Count;

  private:
    // >>> DATA

    // Stored particles
    Stack_Particle d_particles;

    // Number of copies per particle
    Stack_Count d_count;

    // Number of total particles (sum of d_count)
    int d_total;

  public:

    // HOST API

    // Constructor.
    Bank() : d_particles(), d_count(), d_total(0) { /* * */ }

    //! \name std::stack-like operations
    //@{

    //! Is the bank empty?
    bool empty() const
    {
        CHECK(d_particles.size() == d_count.size());
        CHECK(d_particles.empty() ? d_total == 0 : true);
        return d_total == 0;
    }

    //! Return the number of particles left to emit
    size_type size() const { return d_total; }

    //! Just pushing a particle
    void push(const SP_Particle& p, size_type count = 1)
    {
        REQUIRE(p);
        REQUIRE(count > 0);

        // Add a copy of the particle to the stack
        d_particles.push_back(std::make_shared<Particle_t>(*p));
        d_count.push_back(count);

        d_total += count;
    }

    //! Emit the topmost particles from the stack into empty spots in a vector
    void pop( cuda::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles );
    //@}

    //! \name Implementation-specific operations
    //@{

    //! Return the number of unique particles being stored
    size_type num_unique() const { return d_particles.size(); }

    //! Return the number of copies of the next particle
    size_type next_count() const
    {
        return (empty() ? 0 : d_count.back());
    }

    //@}
};

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//

#endif // cuda_mc_Bank_hh

//---------------------------------------------------------------------------//
//              end of cuda_mc/Bank.hh
//---------------------------------------------------------------------------//
